#!/usr/bin/env Rscript
# deg_deseq2.R
#
# Runs a standard DESeq2 differential expression analysis on raw RNA-seq
# count data. Called by pipeline/deg.py when detect.py decides the input
# looks like raw counts (non-negative integers).
#
# Usage:
#   Rscript deg_deseq2.R <counts.csv> <meta.csv> <condition_col> <outdir> [alpha]
#
#   counts.csv  - genes x samples, first column = gene id, header = sample ids
#   meta.csv    - sample id column (must be named "sample") + condition_col
#                 + any extra covariates (used automatically in the design
#                 formula if present)
#   condition_col - name of the column in meta.csv to test (e.g. "diagnosis")
#   outdir      - where results go
#   alpha       - padj cutoff for calling something a DEG (default 0.05)

suppressMessages({
  library(DESeq2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: deg_deseq2.R <counts.csv> <meta.csv> <condition_col> <outdir> [alpha]")
}

counts_path   <- args[1]
meta_path     <- args[2]
condition_col <- args[3]
outdir        <- args[4]
alpha         <- if (length(args) >= 5) as.numeric(args[5]) else 0.05

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_raw <- read.csv(counts_path, row.names = 1, check.names = FALSE)
meta       <- read.csv(meta_path, check.names = FALSE)
rownames(meta) <- meta$sample

# keep only samples present in both files, same order
common <- intersect(colnames(counts_raw), rownames(meta))
if (length(common) < 4) {
  stop("fewer than 4 overlapping samples between counts and metadata - check sample id columns")
}
counts_raw <- counts_raw[, common]
meta <- meta[common, , drop = FALSE]

# round in case of estimated/kallisto-style non-integer counts
counts_mat <- round(as.matrix(counts_raw))
storage.mode(counts_mat) <- "integer"

meta[[condition_col]] <- factor(meta[[condition_col]])

# use any extra covariate columns (age, sex, PMI, batch, RIN, ...) if present
covariates <- setdiff(colnames(meta), c("sample", condition_col))
covariates <- covariates[sapply(meta[covariates], function(x) length(unique(x)) > 1)]

design_formula <- if (length(covariates) > 0) {
  as.formula(paste("~", paste(covariates, collapse = " + "), "+", condition_col))
} else {
  as.formula(paste("~", condition_col))
}
message("DESeq2 design: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = design_formula)
dds <- dds[rowSums(counts(dds) >= 10) >= max(3, floor(ncol(dds) * 0.2)), ]
dds <- DESeq(dds)

res <- results(dds, alpha = alpha)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res_df, file.path(outdir, "deg_results.csv"), row.names = FALSE)

degs <- res_df %>% filter(!is.na(padj), padj < alpha)
write.csv(degs, file.path(outdir, "deg_significant.csv"), row.names = FALSE)

# variance-stabilized matrix - this is what network.R / MEGENA should consume,
# NOT raw counts, since MEGENA correlation assumes roughly continuous/normal data
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(as.data.frame(vst_mat), file.path(outdir, "normalized_expression.csv"))

message(sprintf("DESeq2 done: %d genes tested, %d significant at padj < %.3f",
                 nrow(res_df), nrow(degs), alpha))

# wip note: extend tests for design-formula construction from metadata

# wip note: fix module_deg_pathway_summary.csv output ranking

# wip note: add laptop-memory note on full-genome megena

# wip note: introduce limma normalized-data deg path

# wip note: wire up top-variance-gene + deg restriction for megena

# wip note: extend tests for gse174409 counts matrix loader

# wip note: tighten covariate auto-detection from metadata columns

# wip note: extend tests for deg_results.csv / deg_significant.csv writers

# wip note: cover edge case in module_membership.csv writer

# wip note: refactor batch-correction caveat note

# wip note: correct rscript subprocess runner in utils.py

# wip note: extend tests for meta.csv sample-column validation
