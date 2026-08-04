## =============================================================================
## ModuleViz gallery using bundled example data
##
## Generates:
##   1. longitudinal clustered heatmap
##   2. module line/pattern plots
##   3. side-by-side RNA + ATAC-style heatmap
## =============================================================================

library(ModuleViz)

data("heatmap_example")

out_dir <- file.path("examples", "outputs", "moduleviz_gallery")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
preview_dir <- file.path("inst", "extdata", "results")
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

mat <- heatmap_example$mat
meta <- heatmap_example$meta

## House style for this example dataset: the untreated arm is neutral grey, the
## treated arm is red.  Passing the colours explicitly via `annotation_colors`
## (heatmaps) and `group_colors` (line plots) overrides the palette defaults for
## exactly these levels, so the two conditions read the same way in every panel.
COND_COLS <- c(Control = "grey70", Stimulated = "#C1272D")

## ---- 1. longitudinal heatmap by clustered temporal pattern ------------------
obj <- longitudinal_cluster(
  mat,
  meta,
  sample_col = "SampleID",
  time_col = "TimeHr",
  group_col = "Condition",
  k = 4,
  method = "kmeans",
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 10
)

longitudinal_heatmap(
  obj,
  annotation_cols = c("Condition", "TimeHr"),
  annotation_colors = list(Condition = COND_COLS),
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  top_n_label = 2,
  cluster_palette = "spectral",
  label_fontsize = 11,
  title = "ModuleViz longitudinal heatmap",
  file = file.path(out_dir, "01_longitudinal_clustered_heatmap.pdf"),
  width = 8,
  height = 9
)
longitudinal_heatmap(
  obj,
  annotation_cols = c("Condition", "TimeHr"),
  annotation_colors = list(Condition = COND_COLS),
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  top_n_label = 2,
  cluster_palette = "spectral",
  label_fontsize = 11,
  title = "ModuleViz longitudinal heatmap",
  file = file.path(preview_dir, "01_longitudinal_clustered_heatmap.png"),
  width = 8,
  height = 9
)

## ---- 2. module line/pattern plots ------------------------------------------
p_lines <- pattern_lineplot(
  obj,
  file = file.path(out_dir, "02_module_line_patterns.pdf"),
  group_colors = COND_COLS,
  base_size = 13,
  width = 7,
  height = 6
)
ggplot2::ggsave(
  file.path(preview_dir, "02_module_line_patterns.png"),
  p_lines,
  width = 7,
  height = 6,
  dpi = 300
)

write_memberships(
  obj,
  prefix = file.path(out_dir, "moduleviz_example")
)

## ---- 3. side-by-side RNA + ATAC-style heatmap -------------------------------
## This simulates a paired ATAC layer with two peaks per gene. In a real project,
## use your peak-by-sample ATAC/Cut&Run matrix and row names such as
## "PeakID_GENE".
set.seed(99)
genes <- rownames(mat)
atac <- mat[rep(seq_along(genes), each = 2), , drop = FALSE]
rownames(atac) <- paste0(
  "Peak",
  sprintf("%03d", seq_len(nrow(atac))),
  "_",
  rep(genes, each = 2)
)
atac <- atac + matrix(rnorm(length(atac), sd = 0.35), nrow = nrow(atac))

dual_omics_heatmap(
  primary = obj,
  secondary = atac,
  primary_name = "RNA",
  secondary_name = "ATAC",
  gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = "mean",
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  cluster_palette = "spectral",
  label_fontsize = 11,
  file = file.path(out_dir, "03_rna_atac_side_by_side_heatmap.pdf"),
  width = 11,
  height = 9
)
dual_omics_heatmap(
  primary = obj,
  secondary = atac,
  primary_name = "RNA",
  secondary_name = "ATAC",
  gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = "mean",
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  cluster_palette = "spectral",
  label_fontsize = 11,
  file = file.path(preview_dir, "03_rna_atac_side_by_side_heatmap.png"),
  width = 11,
  height = 9
)

message("Example plots written to: ", normalizePath(out_dir))
message("Preview plots written to: ", normalizePath(preview_dir))
