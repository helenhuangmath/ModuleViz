test_that("example data can be clustered and summarized", {
  data("heatmap_example", package = "HeatmapPlot")

  obj <- longitudinal_cluster(
    heatmap_example$mat,
    heatmap_example$meta,
    sample_col = "SampleID",
    time_col = "TimeHr",
    group_col = "Condition",
    k = 4,
    aggregate_replicates = TRUE,
    stability = TRUE,
    n_resample = 3
  )

  expect_s3_class(obj, "longi")
  expect_equal(nrow(obj$z), nrow(heatmap_example$mat))
  expect_equal(length(levels(obj$cluster)), 4)
  expect_named(obj$membership, c("ID", "Cluster", "MaxAbsZ"))
  expect_named(obj$stability, c("Cluster", "MeanJaccard", "N"))
})

test_that("choose_k returns diagnostics and a suggested k", {
  data("heatmap_example", package = "HeatmapPlot")

  diag <- choose_k(
    zscore_rows(heatmap_example$mat),
    k_range = 2:5,
    method = "kmeans"
  )

  expect_s3_class(diag, "data.table")
  expect_true(attr(diag, "suggested_k") %in% 2:5)
  expect_named(diag, c("k", "tot_withinss", "mean_silhouette", "elbow_k"))
})

test_that("labels resolve from explicit IDs and top features", {
  data("heatmap_example", package = "HeatmapPlot")

  obj <- longitudinal_cluster(
    heatmap_example$mat,
    heatmap_example$meta,
    sample_col = "SampleID",
    time_col = "TimeHr",
    group_col = "Condition",
    k = 4,
    aggregate_replicates = TRUE,
    stability = FALSE
  )

  labels <- resolve_labels(obj, features = "Gene001", top_n = 1)
  expect_s3_class(labels, "data.table")
  expect_true("Gene001" %in% labels$label)
  expect_true(nrow(labels) >= 4)
})
