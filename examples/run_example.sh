#!/usr/bin/env bash
# Full pipeline run against GSE174409 (see data/README.md for how to get it).
set -euo pipefail

python -m pipeline.cli full \
  --expr data/GSE174409_counts.csv \
  --meta data/meta.csv \
  --condition diagnosis \
  --outdir results/oud_dlpfc \
  --top-var 3000 \
  --cores 4 \
  --org-db org.Hs.eg.db

# wip note: document covariate auto-detection from metadata columns

# wip note: fix hypergeometric deg-per-module enrichment

# wip note: clean up pearson correlation network construction
