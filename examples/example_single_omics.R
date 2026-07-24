## =============================================================================
## Example: single-omics longitudinal analysis (e.g. one ATAC or RNA time course)
##
## Mirrors a typical DESeq2 workflow: normalized counts + a DE-stats table are
## used to pick features, then longitudinal_heatmap.R does the clustering and
## the three plots.  Adapt the paths / column names to your project.
## =============================================================================

library(data.table)
library(dplyr)

## Preferred after installing the package:
## install.packages("remotes")
## remotes::install_github("helenhuangmath/Heatmap_Plot")
library(ModuleViz)

## During local development, you can instead source the R files:
# source("../R/longitudinal_heatmap.R")

set.seed(1)

## ---- paths ------------------------------------------------------------------
in_dir  <- "path/to/inputs"
out_dir <- "path/to/Results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(in_dir, "NormRC_DESeq2_filtered.txt")  # feature x sample
stats_file  <- file.path(in_dir, "DESeq2_Wald_pairwise_long_.txt")
meta_file   <- file.path(in_dir, "sample_conditions.txt")       # one row per sample

LFC_CUTOFF  <- 1
PADJ_CUTOFF <- 0.05
K           <- 12
PREFIX      <- file.path(out_dir, sprintf("lfc1_padj05_kmeans%d_%s", K, Sys.Date()))

## ---- 1. load ----------------------------------------------------------------
counts <- fread(counts_file)
setnames(counts, 1, "ID")                       # first column = "PeakID_SYMBOL"

meta <- fread(meta_file)
## meta must have a sample id, a numeric time, and (optionally) a grouping var.
## Here we derive genotype from the Condition string, e.g. "WT_4h" -> "WT".
meta[, Genotype := sub("_.*$", "", Condition)]
## make sure TimeHr is numeric so clusters can be ordered by time
# meta[, TimeHr := as.numeric(TimeHr)]

## ---- 2. choose features to plot (significant genes/peaks) -------------------
stats <- fread(stats_file)
stats[, PeakID_Gene := paste0(PeakID, "_", SYMBOL)]
pass  <- !is.na(stats$padj) & stats$padj < PADJ_CUTOFF &
         !is.na(stats$log2FoldChange) & abs(stats$log2FoldChange) > LFC_CUTOFF
keep_ids <- unique(stats$PeakID_Gene[pass])
counts   <- counts[ID %in% keep_ids]

## optional: a ranking table so the heatmap can auto-label top genes per cluster
rank_tbl <- stats[pass, .(ID = PeakID_Gene,
                          value = -log10(pmax(padj, 1e-300)))]

## ---- 3. cluster into temporal modules --------------------------------------
## (choose_k() first if you want help picking K)
mat_z_input <- as.matrix(counts[, -1]); rownames(mat_z_input) <- counts$ID

obj <- longitudinal_cluster(
  mat_z_input, meta,
  sample_col = "SampleID",
  time_col   = "TimeHr",
  group_col  = "Genotype",     # set NULL if you only have one series
  k = K, method = "kmeans",    # or "hierarchical" / "pam"
  stability = TRUE
)
print(obj)

## ---- 4./5. heatmap with metadata bars + labelled genes ----------------------
longitudinal_heatmap(
  obj,
  annotation_cols = c("Genotype", "TimeHr"),
  top_n_label     = 5,             # auto-label top 5 per cluster ...
  rank_table      = rank_tbl,      # ... ranked by this table
  label_features  = c("Xkr4", "Sox2"),   # ... plus these specific genes
  file   = paste0(PREFIX, "_heatmap.pdf"),
  width  = 10, height = 14
)

## ---- pattern line plot ------------------------------------------------------
pattern_lineplot(obj, file = paste0(PREFIX, "_pattern_lines.pdf"))

## ---- write memberships / centroids / stability ------------------------------
write_memberships(obj, prefix = PREFIX)

## Or do all of the above in one call:
# run_longitudinal(mat_z_input, meta, out_prefix = PREFIX,
#                  sample_col = "SampleID", time_col = "TimeHr",
#                  group_col = "Genotype", k = K,
#                  top_n_label = 5, rank_table = rank_tbl)
