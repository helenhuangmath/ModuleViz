## =============================================================================
## Single-cell example: pseudobulk longitudinal modules
##
## Simulates a small single-cell time course, summarises the cells into
## pseudobulk profiles by sample, cell type, condition, and time point with
## pseudobulk_cells(), and runs the usual ModuleViz workflow on the result.
##
## With real data, replace `sce` by your own SingleCellExperiment or Seurat
## object holding raw counts.
## =============================================================================

library(ModuleViz)

out_dir <- file.path("examples", "outputs", "single_cell_pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
preview_dir <- file.path("inst", "extdata", "results")
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

sce <- simulate_single_cell_example(
  cells_per_sample = 25,
  cell_types = c("T_cell", "Monocyte", "B_cell"),
  seed = 2026
)

## One profile per cell type x biological sample: counts are summed and
## converted to log2 CPM.
pb <- pseudobulk_cells(
  sce,
  time_col = "TimeHr",
  sample_col = "SampleID",
  group_col = "Condition",
  celltype_col = "CellType",
  min_cells = 10
)
pb$Group <- paste(pb$CellType, pb$Condition, sep = "_")

## Condition and CellType are constant within each Group x TimeHr cell, so
## they survive replicate aggregation and stay available as annotation bars.
obj <- longitudinal_cluster(
  pb,
  time_col = "TimeHr",
  group_col = "Group",
  k = 5,
  method = "kmeans",
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 10
)

draw_heatmap <- function(file) {
  longitudinal_heatmap(
    obj,
    ## one bar is enough: the paired palette shades control (light) and
    ## stimulated (dark) within each cell type
    annotation_cols = "Group",
    annotation_palette = "paired",
    top_n_label = 1,
    label_features = c("Gene001", "Gene016", "Gene031"),
    cluster_palette = "spectral",
    file = file,
    width = 7,
    height = 5.5
  )
}
draw_heatmap(file.path(out_dir, "single_cell_pseudobulk_heatmap.pdf"))
draw_heatmap(file.path(preview_dir, "04_single_cell_pseudobulk_heatmap.png"))

p_sc <- pattern_lineplot(
  obj,
  palette = "paired",
  file = file.path(out_dir, "single_cell_pseudobulk_patterns.pdf"),
  width = 5,
  height = 5.5
)
ggplot2::ggsave(
  file.path(preview_dir, "05_single_cell_pseudobulk_patterns.png"),
  p_sc,
  width = 5,
  height = 5.5,
  dpi = 300
)

write_memberships(obj, prefix = file.path(out_dir, "single_cell_pseudobulk"))

message("Single-cell pseudobulk example written to: ", normalizePath(out_dir))
message("Preview plots written to: ", normalizePath(preview_dir))
