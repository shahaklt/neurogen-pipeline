#!/usr/bin/env Rscript
# megena_network.R
#
# Builds a co-expression network with MEGENA (Multiscale Embedded Gene
# Co-expression Network Analysis) and calls modules. Consumes the
# normalized_expression.csv produced by either deg_deseq2.R (vst) or
# deg_limma.R (already normalized), NOT raw counts.
#
# Because MEGENA is O(n^2) on gene-gene correlation, running it on the full
# ~20k-gene expression matrix is slow and mostly noise. By default this
# script restricts the network to the top-N most variable genes plus any
# significant DEGs handed in via --deg-file, which is both faster and more
# biologically defensible for a project this size. Pass --all-genes to
# disable that and run on everything.
#
# Usage:
#   Rscript megena_network.R <normalized_expression.csv> <outdir> \
#       [--deg-file deg_significant.csv] [--top-var 3000] [--cores 4] [--all-genes]

suppressMessages({
  library(MEGENA)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: megena_network.R <normalized_expression.csv> <outdir> [--deg-file f] [--top-var N] [--cores N] [--all-genes]")
}

expr_path <- args[1]
outdir    <- args[2]
rest      <- args[-c(1, 2)]

get_flag <- function(flag, default = NULL) {
  idx <- which(rest == flag)
  if (length(idx) == 0) return(default)
  rest[idx + 1]
}
has_flag <- function(flag) flag %in% rest

deg_file  <- get_flag("--deg-file", NA)
top_var   <- as.numeric(get_flag("--top-var", "3000"))
n_cores   <- as.numeric(get_flag("--cores", "4"))
all_genes <- has_flag("--all-genes")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)

if (!all_genes) {
  keep <- character(0)
  if (!is.na(deg_file) && file.exists(deg_file)) {
    degs <- read.csv(deg_file)
    keep <- intersect(degs$gene, rownames(expr))
    message(sprintf("including %d significant DEGs in the network gene set", length(keep)))
  }
  gene_var <- apply(expr, 1, var)
  var_ranked <- names(sort(gene_var, decreasing = TRUE))
  top <- head(var_ranked, top_var)
  keep <- union(keep, top)
  expr <- expr[rownames(expr) %in% keep, ]
  message(sprintf("network gene set: %d genes (top %d by variance + DEGs)", nrow(expr), top_var))
} else {
  message(sprintf("running on all %d genes (--all-genes) - this will be slow", nrow(expr)))
}

# --- MEGENA standard workflow -------------------------------------------
n.cores <- n_cores
doPar <- n.cores > 1
method <- "pearson"
FDR.cutoff <- 0.05
module.pval <- 0.05
hub.pval <- 0.05
cor.perm <- 10
hub.perm <- 100

ijw <- calculate.correlation(t(expr), doPerm = cor.perm, method = method,
                              FDR.cutoff = FDR.cutoff, n.increment = 100,
                              is.signed = FALSE)

if (doPar) {
  require(doParallel)
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
}

el <- calculate.PFN(ijw[, 1:3], doPar = doPar, num.cores = n.cores)
g <- graph.data.frame(el, directed = FALSE)

MEGENA.output <- do.MEGENA(g,
  mod.pval = module.pval, hub.pval = hub.pval, remove.unsig = TRUE,
  min.size = 10, max.size = vcount(g) / 2,
  doPar = doPar, num.cores = n.cores, n.perm = hub.perm,
  save.output = FALSE)

if (doPar) parallel::stopCluster(cl)

summary.output <- MEGENA.ModuleSummary(MEGENA.output,
  mod.pvalue = module.pval, hub.pvalue = hub.pval,
  min.size = 10, max.size = vcount(g) / 2,
  annot.table = NULL, output.sig = TRUE)

# module membership: gene -> module id, one row per gene
membership <- summary.output$module.table
modules_list <- MEGENA.output$module.output$modules

module_df <- do.call(rbind, lapply(names(modules_list), function(mod_id) {
  data.frame(gene = modules_list[[mod_id]], module = mod_id, stringsAsFactors = FALSE)
}))

write.csv(module_df, file.path(outdir, "module_membership.csv"), row.names = FALSE)
write.csv(el, file.path(outdir, "network_edgelist.csv"), row.names = FALSE)
write.csv(summary.output$module.table, file.path(outdir, "module_summary.csv"), row.names = FALSE)

if (!is.null(MEGENA.output$hub.output$hub.list)) {
  hubs <- do.call(rbind, lapply(names(MEGENA.output$hub.output$hub.list), function(mod_id) {
    hs <- MEGENA.output$hub.output$hub.list[[mod_id]]
    if (length(hs) == 0) return(NULL)
    data.frame(module = mod_id, hub_gene = hs, stringsAsFactors = FALSE)
  }))
  write.csv(hubs, file.path(outdir, "module_hubs.csv"), row.names = FALSE)
}

message(sprintf("MEGENA done: %d genes, %d modules", nrow(module_df),
                 length(unique(module_df$module))))

# wip note: extend tests for planar filtered network step

# wip note: wire up single-cell (pseudo-bulk only) caveat note

# wip note: correct unit tests for detect.py edge cases

# wip note: clarify --all-genes full-transcriptome megena option

# wip note: cover edge case in run_example.sh example script

# wip note: patch brain_region as alternate --condition example

# wip note: fix brain_region as alternate --condition example

# wip note: add vst() variance-stabilizing transform for network input

# wip note: tighten planar filtered network step

# wip note: add brain_region as alternate --condition example

# wip note: simplify cli enrichment subcommand

# wip note: implement module_hubs.csv writer

# wip note: add test for batch-correction caveat note
