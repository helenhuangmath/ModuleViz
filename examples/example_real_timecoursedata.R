## =============================================================================
## Real public dataset example: sorghum drought RNA-seq time course
##
## Data source: Bioconductor package `timecoursedata`.
## This script downloads/loads the public Varoquaux 2019 sorghum leaf RNA-seq
## dataset, keeps a small high-variance subset for speed, and runs the complete
## ModuleViz workflow.
## =============================================================================

library(ModuleViz)

if (!requireNamespace("timecoursedata", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("timecoursedata")
}

out_dir <- file.path("examples", "outputs", "sorghum_leaf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

example <- load_real_timecourse_example(
  dataset = "sorghum_leaf",
  genotype = "BT642",
  min_week = 3,
  top_n_features = 1000
)

message(example$source)
message(example$citation_hint)

diag <- choose_k(
  zscore_rows(example$mat),
  k_range = 2:10,
  method = "kmeans",
  file = file.path(out_dir, "choose_k_elbow.pdf")
)
print(diag)

k <- attr(diag, "suggested_k")

obj <- longitudinal_cluster(
  example$mat,
  example$meta,
  sample_col = example$sample_col,
  time_col = example$time_col,
  group_col = example$group_col,
  k = k,
  method = "kmeans",
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 20
)

print(obj)

longitudinal_heatmap(
  obj,
  annotation_cols = c(example$group_col, example$time_col),
  top_n_label = 3,
  file = file.path(out_dir, "sorghum_leaf_heatmap.pdf"),
  width = 9,
  height = 12
)

pattern_lineplot(
  obj,
  file = file.path(out_dir, "sorghum_leaf_patterns.pdf")
)

write_memberships(
  obj,
  prefix = file.path(out_dir, "sorghum_leaf")
)
