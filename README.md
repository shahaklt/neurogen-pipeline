<p align="center">
  <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg">
  <img alt="Python 3.11" src="https://img.shields.io/badge/python-3.11-blue.svg">
  <img alt="R 4.3" src="https://img.shields.io/badge/R-4.3-276DC3.svg">
  <img alt="DESeq2 / limma / MEGENA" src="https://img.shields.io/badge/DE-DESeq2%20%2F%20limma-9cf.svg">
</p>

<h1 align="center">neurogen-pipeline</h1>

<p align="center">
  Gene expression matrix in, "here's what's different and the biology it points to" out.
</p>

A CLI that takes a gene expression matrix and a sample sheet and runs
differential expression, then builds a MEGENA co-expression network on
top of the DEGs, then runs GO/KEGG enrichment per network module, then
produces one summary table that ties DEG enrichment back to specific
pathways per module. It figures out most of the boring decisions itself
- raw counts vs. already-normalized data, which columns in your
metadata are covariates, what organism you're working in - instead of
making you hand-edit an R script every time the input changes shape.

It was built around one dataset (a postmortem opioid use disorder
RNA-seq study, see below) but nothing about it is hardcoded to that
dataset. Point it at a different counts matrix and a different
`org_db` and it should just work.

---

## Contents

1. [Why this exists](#why-this-exists)
2. [The dataset it was built around](#the-dataset-it-was-built-around)
3. [Setup](#setup)
4. [Run it](#run-it)
5. [How it works](#how-it-works)
6. [Why it's set up this way](#why-its-set-up-this-way)
7. [Output](#output)
8. [What this is and isn't](#what-this-is-and-isnt)
9. [Known limitations](#known-limitations--things-i-havent-gotten-to)
10. [Project layout](#project-layout)
11. [License](#license)

---

## Why this exists

I wanted to reanalyze a public opioid use disorder RNA-seq dataset for
a research project and got tired of rewriting the same
DESeq2 → network → enrichment glue code every time the input format
changed slightly. So instead of one script, this is a small CLI that
looks at whatever expression matrix you give it and picks the right
analysis path itself.

## The dataset it was built around

[GSE174409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409)
- postmortem human dorsolateral prefrontal cortex (DLPFC) and nucleus
accumbens (NAc) RNA-seq, opioid use disorder (n=20) vs. matched
comparison subjects (n=20), from:

> Seney ML, Kim SM, Glausier JR, et al. *Transcriptional Alterations in
> Dorsolateral Prefrontal Cortex and Nucleus Accumbens Implicate
> Neuroinflammation and Synaptic Remodeling in Opioid Use Disorder.*
> Biological Psychiatry, 2021.

Picked this one because it's public, has a real clinical diagnosis
column to test against, and the neuroinflammation/synaptic remodeling
angle in the original paper gives something concrete to check the
pipeline against - if the modules that come out aren't at least
touching immune/synaptic GO terms, something's wrong upstream.

Nothing from this dataset is checked into the repo (see
`data/README.md`) - raw counts and metadata for real subjects shouldn't
live in git history even for a de-identified public study.

---

## Setup

```bash
conda env create -f environment.yml
conda activate neurogen-pipeline

# MEGENA isn't on conda-forge/bioconda, install it separately:
R -e 'install.packages("MEGENA", repos = "https://cloud.r-project.org")'
```

---

## Run it

Whole thing, one command:

```bash
python -m pipeline.cli full \
  --expr data/GSE174409_counts.csv \
  --meta data/meta.csv \
  --condition diagnosis \
  --outdir results/oud_dlpfc
```

Or stage by stage - useful when you're iterating on the network step
and don't want to rerun DESeq2 every time:

```bash
python -m pipeline.cli deg       --expr counts.csv --meta meta.csv --condition diagnosis --outdir results/deg
python -m pipeline.cli network   --expr results/deg/normalized_expression.csv --deg-file results/deg/deg_significant.csv --outdir results/network
python -m pipeline.cli enrichment --modules results/network/module_membership.csv --outdir results/enrichment
```

`meta.csv` needs a `sample` column matching your expression matrix's
column headers, plus whatever condition column you want to test
(`--condition diagnosis`, `--condition brain_region`, whatever). Any
other columns (age, sex, PMI, RIN, batch...) get folded into the model
formula automatically as covariates - that's the "adaptive" part, see
below.

Different organism, e.g. a mouse model instead of human postmortem
tissue - swap the `OrgDb`:

```bash
python -m pipeline.cli enrichment --modules module_membership.csv \
  --org-db org.Mm.eg.db --outdir results/enrichment
```

---

## How it works

```
counts.csv + meta.csv
        │
        ▼
 ┌───────────────┐   detect.py looks at the values: integer counts -> DESeq2,
 │   DEG step     │   anything else (microarray, logCPM, ...) -> limma.
 │ deseq2/limma   │   Extra metadata columns become covariates automatically.
 └───────┬───────┘
         │ significant DEGs + variance-stabilized expression matrix
         ▼
 ┌───────────────┐   Pearson correlation network, planar filtered network,
 │  MEGENA step   │   multiscale module detection. Restricted to top-variance
 │                │   genes + DEGs by default (full-transcriptome MEGENA is
 └───────┬───────┘   slow and mostly noise for a dataset this size).
         │ module membership + hub genes
         ▼
 ┌───────────────┐   clusterProfiler GO (BP) + KEGG over-representation,
 │  GO/KEGG step  │   run separately per module.
 └───────┬───────┘
         │
         ▼
 module_deg_pathway_summary.csv
   (hypergeometric DEG enrichment per module + its top GO terms,
    sorted so the modules most worth reading about are at the top)
```

`pipeline/` is the orchestration layer (Python); it calls out to the
actual analysis code in `R/` (DESeq2, limma, MEGENA, clusterProfiler)
via `Rscript`.

---

## Why it's set up this way

The "adaptive" part is mostly in `pipeline/detect.py`. It's not fancy -
it checks whether the values in your matrix are non-negative integers
(raw counts → DESeq2) or not (assume it's already normalized → limma),
and it pulls in whatever extra columns are in your metadata file as
covariates in the model formula. That's enough to reuse this on a
microarray dataset or someone else's RNA-seq counts without touching
the R scripts, which is the whole point - I didn't want to hand-edit a
design formula every time I pointed this at a different GEO series.

**MEGENA gets the DESeq2 `vst()`-transformed matrix, not raw counts.**
Correlation on raw counts is dominated by library size effects - a
mistake I made the first time I ran this, before spending an afternoon
figuring out why every module looked like it was just "high expression
genes."

**By default MEGENA only runs on the top 3000 most-variable genes plus
whatever came out of the DEG step, not everything.** The full-genome
correlation matrix on ~20k genes made my laptop fall over the first
time and mostly finds noise modules on a dataset with 40 samples
anyway. Pass `--all-genes` if you actually want that.

**GO/KEGG enrichment runs per module, not once over the whole DEG
list.** A flat enrichment table over all DEGs tells you the dataset
involves immune and synaptic biology; it doesn't tell you which
co-expression module is driving which piece of that. Running it per
module and then correlating module membership against the DEG list is
what makes `module_deg_pathway_summary.csv` the file worth actually
reading.

---

## Output

- `deg/deg_results.csv`, `deg/deg_significant.csv` - full and filtered DE results
- `deg/normalized_expression.csv` - vst/normalized matrix (feeds the network step)
- `network/module_membership.csv` - gene → module assignment
- `network/module_hubs.csv` - hub genes per module
- `network/network_edgelist.csv` - full planar filtered network, importable into Cytoscape
- `enrichment/go_enrichment_by_module.csv`, `enrichment/kegg_enrichment_by_module.csv`
- `module_deg_pathway_summary.csv` - the one to actually read first

---

## What this is and isn't

Read this before pointing it at a dataset you're publishing on.

**This is a real DESeq2/limma → MEGENA → clusterProfiler pipeline.**
Every stage calls the actual, standard Bioconductor/CRAN packages you'd
use if you wrote this by hand - nothing here reimplements DE testing,
network construction, or enrichment statistics. The Python layer is
orchestration and format detection, not analysis.

**The "adaptive" detection is a couple of heuristics, not a model.**
Counts-vs-normalized detection is "are all values non-negative
integers." Covariate handling is "every metadata column that isn't
`sample` or your condition column goes in the formula." That's enough
to reuse across datasets with the same shape of problem; it is not a
substitute for looking at your design matrix before you trust the
output.

**It was validated against exactly one dataset (GSE174409).** The
detection logic has unit tests (`tests/test_detect.py`); the DE →
network → enrichment path end to end has been run against this one
study, not benchmarked across many. Treat a first run on a new dataset
as something to sanity-check, not something to trust blind.

**It doesn't do batch correction, multi-group contrasts, or
single-cell.** See [Known limitations](#known-limitations--things-i-havent-gotten-to)
for specifics - none of these are silently handled for you.

---

## Known limitations / things I haven't gotten to

- `config/example_config.yaml` isn't wired up yet - everything's still
  CLI flags. Wanted a `--config file.yaml` option so I stop retyping
  the same flags but haven't done it.
- No multi-group (>2 level) contrast handling beyond whatever
  DESeq2/limma do by default with the reference level - if your
  condition column has 3+ groups you'll want to specify contrasts
  yourself.
- Batch correction (ComBat, etc.) isn't part of the pipeline. If your
  dataset needs it, run it before this and feed in the corrected
  matrix.
- No single-cell support. This assumes bulk (pseudo-bulked, if it's
  originally single-cell).
- Only tested end to end on GSE174409. The DEG step has unit tests for
  the detection logic; I haven't set up fixtures for the R side
  (MEGENA in particular is slow enough that a CI run isn't really
  practical on free-tier runners).

---

## Project layout

```
pipeline/          adaptive orchestration layer (python)
  detect.py         figures out what kind of data you gave it
  deg.py            calls the right R script for DE
  network.py        calls MEGENA
  enrichment.py     calls GO/KEGG enrichment
  cli.py            argparse entrypoint
  utils.py          Rscript runner + module/DEG/pathway correlation
R/                 the actual analysis code (DESeq2, limma, MEGENA, clusterProfiler)
data/              not checked in, see data/README.md
tests/             unit tests for the python detection logic
examples/          example run script
config/            documents the CLI's knobs; not wired up yet
```

---

## License

MIT. Use it, change it, ship it, cite it, fork it and never mention this
repo again - no permission needed and no attribution expected beyond
keeping the license file with the code.
