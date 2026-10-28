## =============================================================================
## ModuleViz gallery using bundled example data
##
## Generates:
##   1. longitudinal clustered heatmap
##   2. module line/pattern plots
##   3. side-by-side RNA + ATAC-style heatmap (dual_heatmap)
##   4. RNA + ATAC + H3K27ac + H3K27me3 heatmap (multi_heatmap)
##   5. module trajectories in every layer (multi_lineplot)
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
## use your peak-by-sample ATAC/CUT&RUN matrix and row names such as
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

dual_heatmap(
  primary = obj,
  secondary = atac,
  primary_name = "RNA",
  secondary_name = "ATAC",
  gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = "mean",
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  cluster_palette = "spectral",
  file = file.path(out_dir, "03_rna_atac_side_by_side_heatmap.pdf"),
  width = 11,
  height = 9
)
dual_heatmap(
  primary = obj,
  secondary = atac,
  primary_name = "RNA",
  secondary_name = "ATAC",
  gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = "mean",
  label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
  cluster_palette = "spectral",
  file = file.path(preview_dir, "03_rna_atac_side_by_side_heatmap.png"),
  width = 11,
  height = 9
)

## ---- 4. several layers side by side: RNA + ATAC + H3K27ac + H3K27me3 ------
## Two simulated histone layers: an activating mark (H3K27ac) that follows
## expression and a repressive mark (H3K27me3) that mirrors it.  Each layer is
## z-scored within itself and collapsed to one row per gene, then drawn as its
## own panel in the RNA module order.
as_regions <- function(m, sign = 1, sd = 0.5) {
  out <- sign * m + matrix(rnorm(length(m), sd = sd), nrow = nrow(m))
  rownames(out) <- paste0("region", seq_len(nrow(m)), "_", rownames(m))
  out
}
layers <- list(RNA      = obj,
               ATAC     = atac,
               H3K27ac  = as_regions(mat),
               H3K27me3 = as_regions(mat, sign = -1))

multi <- multi_cluster(layers, meta = meta, sample_col = "SampleID",
                       time_col = "TimeHr", group_col = "Condition")
for (f in c(file.path(out_dir, "06_multi_heatmap.pdf"),
            file.path(preview_dir, "06_multi_heatmap.png"))) {
  multi_heatmap(
    multi,
    annotation_colors = list(Condition = COND_COLS),
    label_features = c("Gene001", "Gene016", "Gene031", "Gene046"),
    cluster_palette = "spectral",
    show_column_names = FALSE,
    file = f,
    height = 9
  )
}

## ---- 5. module trajectories in every layer ----------------------------------
for (f in c(file.path(out_dir, "07_multi_lineplot.pdf"),
            file.path(preview_dir, "07_multi_lineplot.png"))) {
  multi_lineplot(multi, base_size = 13, file = f, width = 8, height = 8)
}

message("Example plots written to: ", normalizePath(out_dir))
message("Preview plots written to: ", normalizePath(preview_dir))
