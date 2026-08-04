## =============================================================================
## Example: two data sets side by side (RNA + ATAC, or RNA + Cut&Run)
##
## The primary layer (RNA) drives the clustering and row order; the secondary
## layer (ATAC / Cut&Run) is aligned per gene and shown alongside so you can see
## expression and chromatin signal for the same genes together.  This produces a
## SEPARATE new figure, independent of the single-omics heatmaps.
## =============================================================================

library(data.table)

## Preferred after installing the package:
## install.packages("remotes")
## remotes::install_github("helenhuangmath/ModuleViz")
library(ModuleViz)

## During local development, you can instead source the R files:
# source("../R/longitudinal_heatmap.R")
set.seed(1)

out_dir <- "path/to/Results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- load both layers -------------------------------------------------------
## RNA: one row per gene (rownames = gene symbol), columns = samples.
rna_counts <- fread("path/to/RNA_NormRC_filtered.txt")
rna <- as.matrix(rna_counts[, -1]); rownames(rna) <- rna_counts[[1]]  # gene symbols

## ATAC / Cut&Run: rows = peaks named "PeakID_SYMBOL", columns = samples.
atac_counts <- fread("path/to/ATAC_NormRC_filtered.txt")
atac <- as.matrix(atac_counts[, -1]); rownames(atac) <- atac_counts[[1]]

## Sample metadata (shared or per-layer). Must map columns -> time (+ group).
meta <- fread("path/to/sample_conditions.txt")
meta[, Genotype := sub("_.*$", "", Condition)]

## ---- Option A: cluster RNA first, then attach ATAC --------------------------
rna_obj <- longitudinal_cluster(
  rna, meta,
  sample_col = "SampleID", time_col = "TimeHr", group_col = "Genotype",
  k = 10, method = "kmeans", stability = TRUE
)

dual_omics_heatmap(
  primary        = rna_obj,
  secondary      = atac,
  primary_name   = "RNA",
  secondary_name = "ATAC",
  ## how to turn a feature ID into a gene symbol:
  gene_of_primary   = function(id) id,                 # RNA rownames are symbols
  gene_of_secondary = function(id) sub("^.*_", "", id),# "PeakID_SYMBOL" -> "SYMBOL"
  aggregate      = "mean",     # collapse multiple peaks per gene by mean (or "top")
  label_features = c("Sox2", "Pax6", "Gfap"),
  file   = file.path(out_dir, sprintf("RNA_ATAC_sidebyside_%s.pdf", Sys.Date())),
  width  = 12, height = 14
)

## ---- Option B: pass raw matrices and let the function cluster the primary ----
## (clustering args are forwarded to longitudinal_cluster via ...)
# dual_omics_heatmap(
#   primary = rna, secondary = atac,
#   primary_name = "RNA", secondary_name = "CutRun",
#   meta = meta, sample_col = "SampleID", time_col = "TimeHr",
#   group_col = "Genotype", k = 10,
#   gene_of_secondary = function(id) sub("^.*_", "", id),
#   aggregate = "top",
#   file = file.path(out_dir, "RNA_CutRun_sidebyside.pdf")
# )
