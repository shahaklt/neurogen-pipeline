#!/usr/bin/env Rscript
# deg_limma.R
#
# Differential expression for data that's already normalized/continuous -
# microarray intensities, log2 CPM, proteomics, whatever isn't raw integer
# counts. Called by pipeline/deg.py when detect.py flags the input as
# non-count data.
#
# Usage:
#   Rscript deg_limma.R <expr.csv> <meta.csv> <condition_col> <outdir> [alpha]

suppressMessages({
  library(limma)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: deg_limma.R <expr.csv> <meta.csv> <condition_col> <outdir> [alpha]")
}

expr_path     <- args[1]
meta_path     <- args[2]
condition_col <- args[3]
outdir        <- args[4]
alpha         <- if (length(args) >= 5) as.numeric(args[5]) else 0.05

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)
rownames(meta) <- meta$sample

common <- intersect(colnames(expr), rownames(meta))
if (length(common) < 4) {
  stop("fewer than 4 overlapping samples between expression matrix and metadata")
}
expr <- as.matrix(expr[, common])
meta <- meta[common, , drop = FALSE]
meta[[condition_col]] <- factor(meta[[condition_col]])

covariates <- setdiff(colnames(meta), c("sample", condition_col))
covariates <- covariates[sapply(meta[covariates], function(x) length(unique(x)) > 1)]

design_terms <- c(covariates, condition_col)
design <- model.matrix(as.formula(paste("~", paste(design_terms, collapse = " + "))), data = meta)
message("limma design columns: ", paste(colnames(design), collapse = ", "))

fit <- lmFit(expr, design)
fit <- eBayes(fit)

coef_name <- grep(condition_col, colnames(design), value = TRUE)[1]
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("gene")

write.csv(tt, file.path(outdir, "deg_results.csv"), row.names = FALSE)

degs <- tt %>% filter(adj.P.Val < alpha)
write.csv(degs, file.path(outdir, "deg_significant.csv"), row.names = FALSE)

# expr matrix is already normalized/continuous - pass it straight through
# for the network step
write.csv(as.data.frame(expr), file.path(outdir, "normalized_expression.csv"))

message(sprintf("limma done: %d features tested, %d significant at adj.P < %.3f",
                 nrow(tt), nrow(degs), alpha))

# wip note: tighten hypergeometric deg-per-module enrichment
