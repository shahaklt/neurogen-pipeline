# neurogen-pipeline

> **Note on this repo's commit history:** the commits in this repo were
> generated in one batch for a classroom demonstration of git
> commit-date spoofing — every author/committer date was set explicitly
> via `GIT_AUTHOR_DATE`/`GIT_COMMITTER_DATE`, none of them reflects when
> the corresponding change was actually written. The code itself is a
> real sample project; the *history* is a teaching prop. See
> [shahaklt/git-commit-date-spoofing-demo](https://github.com/shahaklt/git-commit-date-spoofing-demo)
> for the mechanism behind it.

A pipeline for going from a gene expression matrix to "here's what's
actually different and here's the biology it points to" - differential
expression, then a MEGENA co-expression network built on top of the
DEGs, then GO/KEGG enrichment per network module, then a summary table
that ties DEG enrichment back to specific pathways per module.

I started this because I wanted to reanalyze a public opioid use disorder
RNA-seq dataset for a research project and got tired of rewriting the
same DESeq2 -> network -> enrichment glue code every time the input
format changed slightly. So instead of one script, this is a small CLI
that looks at whatever expression matrix you give it and picks the right
analysis path itself (raw counts vs already-normalized data, what
covariates exist, etc). It's built around one dataset but isn't hardcoded
to it.

## The dataset

[GSE174409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409) -
postmortem human dorsolateral prefrontal cortex (DLPFC) and nucleus
accumbens (NAc) RNA-seq, opioid use disorder (n=20) vs matched comparison
subjects (n=20), from:

> Seney ML, Kim SM, Glausier JR, et al. *Transcriptional Alterations in
> Dorsolateral Prefrontal Cortex and Nucleus Accumbens Implicate
> Neuroinflammation and Synaptic Remodeling in Opioid Use Disorder.*
> Biological Psychiatry, 2021.

I picked this one because it's public, has a real clinical
diagnosis column to test against, and the neuroinflammation/synaptic
remodeling angle in the original paper gives something concrete to check
the pipeline against - if the modules that come out aren't at least
touching immune/synaptic GO terms, something's wrong upstream. See
`data/README.md` for how to pull the counts matrix.

## What it does

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

## Why it's set up this way

The "adaptive" part is mostly in `pipeline/detect.py`. It's not fancy -
it checks whether the values in your matrix are non-negative integers
(raw counts -> DESeq2) or not (assume it's already normalized -> limma),
and it pulls in whatever extra columns are in your metadata file as
covariates in the model formula. That's enough to reuse this on a
microarray dataset or someone else's RNA-seq counts without touching the
R scripts, which is the whole point - I didn't want to hand-edit a design
formula every time I pointed this at a different GEO series.

MEGENA gets the DESeq2 `vst()`-transformed matrix, not raw counts -
correlation on raw counts is dominated by library size effects, and
that's a mistake I made the first time I ran this before duct hunting
down why every module looked like it was just "high expression genes."

By default MEGENA only runs on the top 3000 most-variable genes plus
whatever came out of the DEG step, not everything - the full-genome
correlation matrix on ~20k genes made my laptop fall over the first time
and mostly finds noise modules on a dataset with 40 samples anyway. Pass
`--all-genes` if you actually want that.

## Setup

```bash
conda env create -f environment.yml
conda activate neurogen-pipeline
# MEGENA isn't on conda-forge, install separately:
R -e 'install.packages("MEGENA", repos = "https://cloud.r-project.org")'
```

## Usage

Run the whole thing:

```bash
python -m pipeline.cli full \
  --expr data/GSE174409_counts.csv \
  --meta data/meta.csv \
  --condition diagnosis \
  --outdir results/oud_dlpfc
```

Or run stages independently (useful when you're iterating on the network
step and don't want to rerun DESeq2 every time):

```bash
python -m pipeline.cli deg --expr counts.csv --meta meta.csv --condition diagnosis --outdir results/deg
python -m pipeline.cli network --expr results/deg/normalized_expression.csv --deg-file results/deg/deg_significant.csv --outdir results/network
python -m pipeline.cli enrichment --modules results/network/module_membership.csv --outdir results/enrichment
```

`meta.csv` needs a `sample` column matching your expression matrix's
column headers and whatever condition column you want to test
(`--condition diagnosis`, `--condition brain_region`, whatever). Any
other columns (age, sex, PMI, RIN, batch...) get folded into the model
formula automatically as covariates.

To point this at a different organism, e.g. a mouse model instead of
human postmortem tissue:

```bash
python -m pipeline.cli enrichment --modules module_membership.csv \
  --org-db org.Mm.eg.db --outdir results/enrichment
```

## Output

- `deg/deg_results.csv`, `deg/deg_significant.csv` - full and filtered DE results
- `deg/normalized_expression.csv` - vst/normalized matrix (this is what feeds the network step)
- `network/module_membership.csv` - gene -> module assignment
- `network/module_hubs.csv` - hub genes per module
- `network/network_edgelist.csv` - full planar filtered network, importable into Cytoscape
- `enrichment/go_enrichment_by_module.csv`, `enrichment/kegg_enrichment_by_module.csv`
- `module_deg_pathway_summary.csv` - the one to actually read first

## Known limitations / things I haven't gotten to

- `config/example_config.yaml` isn't actually wired up yet - everything's
  still CLI flags. Wanted a `--config file.yaml` option so I stop
  retyping the same flags but haven't done it.
- No multi-group (>2 level) contrast handling beyond whatever DESeq2/limma
  do by default with the reference level - if your condition column has
  3+ groups you'll want to specify contrasts yourself.
- Batch correction (ComBat, etc.) isn't part of the pipeline. If your
  dataset needs it, run it before this and feed in the corrected matrix.
- No single-cell support. This assumes bulk (pseudo-bulked, if it's
  originally single-cell).
- Only tested end to end on GSE174409. The DEG step has unit tests for
  the detection logic; I haven't set up fixtures for the R side (MEGENA
  in particular is slow enough that a CI run isn't really practical on
  free-tier runners).

## Repo layout

```
pipeline/          adaptive orchestration layer (python)
  detect.py         figures out what kind of data you gave it
  deg.py            calls the right R script for DE
  network.py        calls MEGENA
  enrichment.py     calls GO/KEGG enrichment
  cli.py            argparse entrypoint
  utils.py          Rscript runner + module/DEG/pathway correlation
R/                 the actual analysis code (DESeq2, limma, MEGENA, clusterProfiler)
data/               not checked in, see data/README.md
tests/              unit tests for the python detection logic
examples/           example run script
```
