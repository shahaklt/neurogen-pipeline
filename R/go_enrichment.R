#!/usr/bin/env Rscript
# go_enrichment.R
#
# Runs GO (and KEGG) over-representation on every MEGENA module, using
# clusterProfiler + org.Hs.eg.db (swap in another org.*.eg.db for
# non-human data - see --org-db).
#
# Usage:
#   Rscript go_enrichment.R <module_membership.csv> <outdir> \
#       [--org-db org.Hs.eg.db] [--key-type SYMBOL] [--qval 0.05] [--universe expr.csv]

suppressMessages({
  library(clusterProfiler)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: go_enrichment.R <module_membership.csv> <outdir> [--org-db pkg] [--key-type SYMBOL] [--qval 0.05] [--universe expr.csv]")
}

module_path <- args[1]
outdir      <- args[2]
rest        <- args[-c(1, 2)]

get_flag <- function(flag, default = NULL) {
  idx <- which(rest == flag)
  if (length(idx) == 0) return(default)
  rest[idx + 1]
}

org_pkg   <- get_flag("--org-db", "org.Hs.eg.db")
key_type  <- get_flag("--key-type", "SYMBOL")
qval_cut  <- as.numeric(get_flag("--qval", "0.05"))
universe_path <- get_flag("--universe", NA)

suppressMessages(library(org_pkg, character.only = TRUE))
org_db <- get(org_pkg)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

modules <- read.csv(module_path)
universe <- if (!is.na(universe_path) && file.exists(universe_path)) {
  rownames(read.csv(universe_path, row.names = 1, check.names = FALSE))
} else {
  NULL
}

mod_ids <- unique(modules$module)
message(sprintf("running GO + KEGG enrichment on %d modules", length(mod_ids)))

go_all <- list()
kegg_all <- list()

for (mod_id in mod_ids) {
  genes <- modules$gene[modules$module == mod_id]
  if (length(genes) < 5) next

  go_res <- tryCatch(
    enrichGO(gene = genes, OrgDb = org_db, keyType = key_type, ont = "BP",
             universe = universe, pAdjustMethod = "BH",
             qvalueCutoff = qval_cut, readable = FALSE),
    error = function(e) NULL)

  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    df <- as.data.frame(go_res)
    df$module <- mod_id
    go_all[[mod_id]] <- df
  }

  # KEGG needs entrez ids
  entrez <- tryCatch(
    suppressMessages(bitr(genes, fromType = key_type, toType = "ENTREZID", OrgDb = org_db)$ENTREZID),
    error = function(e) character(0))

  if (length(entrez) >= 5) {
    kegg_res <- tryCatch(
      enrichKEGG(gene = entrez, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = qval_cut),
      error = function(e) NULL)
    if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
      df <- as.data.frame(kegg_res)
      df$module <- mod_id
      kegg_all[[mod_id]] <- df
    }
  }
}

go_combined <- if (length(go_all) > 0) bind_rows(go_all) else data.frame()
kegg_combined <- if (length(kegg_all) > 0) bind_rows(kegg_all) else data.frame()

write.csv(go_combined, file.path(outdir, "go_enrichment_by_module.csv"), row.names = FALSE)
write.csv(kegg_combined, file.path(outdir, "kegg_enrichment_by_module.csv"), row.names = FALSE)

message(sprintf("GO enrichment done: %d modules with >=1 significant BP term, %d with KEGG hits",
                 length(go_all), length(kegg_all)))

# wip note: extend tests for --org-db flag for non-human organisms

# wip note: note known-limitations section in readme

# wip note: add docstring for run_example.sh example script

# wip note: add docstring for library-size caveat on raw-count correlation

# wip note: fix --org-db flag for non-human organisms

# wip note: tighten normalized_expression.csv writer

# wip note: refactor run_example.sh example script

# wip note: add run_example.sh example script

# wip note: implement network_edgelist.csv cytoscape export
