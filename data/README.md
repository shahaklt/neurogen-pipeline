# data/

Nothing is checked into this folder on purpose (see `.gitignore`) - raw
counts and sample metadata for real subjects shouldn't live in git history
even if the study is technically public/de-identified.

## Dataset this project was built around

**GSE174409** - RNA-seq from postmortem human dorsolateral prefrontal
cortex (DLPFC) and nucleus accumbens (NAc), opioid use disorder (OUD, n=20)
vs matched unaffected comparison subjects (n=20).

Paper: Seney ML, Kim SM, Glausier JR, et al. *Transcriptional Alterations
in Dorsolateral Prefrontal Cortex and Nucleus Accumbens Implicate
Neuroinflammation and Synaptic Remodeling in Opioid Use Disorder.*
Biological Psychiatry, 2021.

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174409

### Getting it

```bash
# raw counts + series metadata
wget -O data/GSE174409_counts.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174409/suppl/GSE174409_raw_counts.txt.gz"
gunzip data/GSE174409_counts.txt.gz
```

Then build `meta.csv` from the GEO sample sheet (`sample`, `diagnosis`,
`brain_region`, plus whatever covariates you want to control for - age,
sex, PMI, RIN all show up in the series matrix). I mostly did this by hand
in a spreadsheet the first time through, there's no script for it because
every GEO series formats its sample table slightly differently.

Expected shapes for the pipeline:
- `counts.csv`: genes (rows) x samples (cols), first column = gene symbol
- `meta.csv`: one row per sample, must have a `sample` column matching the
  counts column headers exactly

That's it - the pipeline figures out counts-vs-normalized and builds the
model formula from whatever other columns are in meta.csv.

# wip note: introduce planar filtered network step

# wip note: simplify module_hubs.csv writer
