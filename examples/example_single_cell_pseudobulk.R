## =============================================================================
## Single-cell-style example: pseudobulk longitudinal modules
##
## This example simulates cell-level expression from the bundled matrix, then
## pseudobulks cells by cell type, condition, time point, and replicate before
## running the same ModuleViz workflow.
## =============================================================================

library(ModuleViz)

data("heatmap_example")

out_dir <- file.path("examples", "outputs", "single_cell_pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bulk_mat <- heatmap_example$mat
bulk_meta <- heatmap_example$meta
genes <- rownames(bulk_mat)
cell_types <- c("T_cell", "Monocyte", "B_cell")
cells_per_sample_type <- 25

set.seed(2026)
cell_meta <- do.call(rbind, lapply(seq_len(nrow(bulk_meta)), function(i) {
  do.call(rbind, lapply(cell_types, function(ct) {
    data.frame(
      CellID = sprintf("%s_%s_Cell%03d", bulk_meta$SampleID[i], ct,
                       seq_len(cells_per_sample_type)),
      SampleID = bulk_meta$SampleID[i],
      Condition = bulk_meta$Condition[i],
      TimeHr = bulk_meta$TimeHr[i],
      Replicate = bulk_meta$Replicate[i],
      CellType = ct,
      stringsAsFactors = FALSE
    )
  }))
}))

cell_mat <- matrix(NA_real_, nrow = nrow(bulk_mat), ncol = nrow(cell_meta),
                   dimnames = list(genes, cell_meta$CellID))
cell_type_effect <- c(T_cell = 0.35, Monocyte = -0.15, B_cell = 0.05)

for (j in seq_len(nrow(cell_meta))) {
  sample_id <- cell_meta$SampleID[j]
  ct <- cell_meta$CellType[j]
  mu <- bulk_mat[, sample_id] + cell_type_effect[[ct]]
  cell_mat[, j] <- mu + rnorm(length(mu), sd = 0.55)
}

pseudobulk_group <- paste(cell_meta$CellType, cell_meta$Condition,
                          cell_meta$TimeHr, cell_meta$Replicate, sep = "|")
pseudobulk_ids <- unique(pseudobulk_group)

pseudobulk_mat <- sapply(pseudobulk_ids, function(id) {
  rowMeans(cell_mat[, pseudobulk_group == id, drop = FALSE])
})
pseudobulk_mat <- as.matrix(pseudobulk_mat)
rownames(pseudobulk_mat) <- genes
colnames(pseudobulk_mat) <- paste0("PB", seq_along(pseudobulk_ids))

parts <- strsplit(pseudobulk_ids, "\\|")
pseudobulk_meta <- data.frame(
  SampleID = colnames(pseudobulk_mat),
  CellType = vapply(parts, `[`, character(1), 1),
  Condition = vapply(parts, `[`, character(1), 2),
  TimeHr = as.numeric(vapply(parts, `[`, character(1), 3)),
  Replicate = as.integer(vapply(parts, `[`, character(1), 4)),
  stringsAsFactors = FALSE
)
pseudobulk_meta$Group <- paste(pseudobulk_meta$CellType,
                               pseudobulk_meta$Condition, sep = "_")

obj <- longitudinal_cluster(
  pseudobulk_mat,
  pseudobulk_meta,
  sample_col = "SampleID",
  time_col = "TimeHr",
  group_col = "Group",
  k = 5,
  method = "kmeans",
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 10
)

longitudinal_heatmap(
  obj,
  annotation_cols = c("Group", "TimeHr"),
  top_n_label = 2,
  label_features = c("Gene001", "Gene016", "Gene031"),
  label_fontsize = 11,
  cluster_palette = "spectral",
  title = "Single-cell pseudobulk temporal modules",
  file = file.path(out_dir, "single_cell_pseudobulk_heatmap.pdf"),
  width = 11,
  height = 10
)

pattern_lineplot(
  obj,
  palette = "science",
  base_size = 12,
  file = file.path(out_dir, "single_cell_pseudobulk_patterns.pdf"),
  width = 9,
  height = 7
)

write_memberships(obj, prefix = file.path(out_dir, "single_cell_pseudobulk"))

message("Single-cell pseudobulk example written to: ", normalizePath(out_dir))
