test_that("example data can be clustered and summarized", {
  data("heatmap_example", package = "ModuleViz")

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
  data("heatmap_example", package = "ModuleViz")

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
  data("heatmap_example", package = "ModuleViz")

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

test_that("pattern line plot accepts journal palette choices", {
  data("heatmap_example", package = "ModuleViz")

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

  p_nature <- pattern_lineplot(obj, palette = "nature")
  p_science <- pattern_lineplot(obj, palette = "science")

  expect_s3_class(p_nature, "ggplot")
  expect_s3_class(p_science, "ggplot")
})

multi_layers <- function() {
  data("heatmap_example", package = "ModuleViz")
  mat <- heatmap_example$mat
  as_peaks <- function(m, per, drop = integer(0)) {
    g <- rownames(m)[setdiff(seq_len(nrow(m)), drop)]
    out <- m[rep(g, each = per), ]
    rownames(out) <- paste0("peak", seq_len(nrow(out)), "_", rep(g, each = per))
    out
  }
  list(RNA = mat, ATAC = as_peaks(mat, 2, drop = 1:3),
       H3K27ac = as_peaks(mat, 3), H3K27me3 = as_peaks(-mat, 1, drop = 4))
}

test_that("multi_cluster aligns every layer to the clustering layer", {
  data("heatmap_example", package = "ModuleViz")
  layers <- multi_layers()

  x <- multi_cluster(layers, meta = heatmap_example$meta,
                           sample_col = "SampleID", time_col = "TimeHr",
                           group_col = "Condition", k = 4, stability = FALSE)

  expect_s3_class(x, "multi_longi")
  expect_named(x$layers, names(layers))
  ## genes 1-4 are missing from ATAC or H3K27me3, so the intersection drops them
  expect_equal(length(x$genes), nrow(layers$RNA) - 4)
  for (L in x$layers) expect_identical(rownames(L$z), x$membership$ID)
  expect_false(anyNA(x$layers$ATAC$z))
  expect_named(x$membership, c("ID", "Gene", "Cluster", "MaxAbsZ"))

  ## gene_set = "reference" keeps every RNA gene and leaves gaps as NA
  r <- multi_cluster(layers, meta = heatmap_example$meta,
                           sample_col = "SampleID", time_col = "TimeHr",
                           k = 4, stability = FALSE, gene_set = "reference")
  expect_equal(length(r$genes), nrow(layers$RNA))
  expect_equal(sum(rowSums(is.na(r$layers$ATAC$z)) > 0), 3)
})

test_that("joint clustering and the multi-omics plots run", {
  data("heatmap_example", package = "ModuleViz")

  f <- tempfile(fileext = ".pdf")
  x <- multi_heatmap(multi_layers(), meta = heatmap_example$meta,
                           sample_col = "SampleID", time_col = "TimeHr",
                           group_col = "Condition", cluster_by = "joint",
                           k = 4, stability = FALSE,
                           label_features = "Gene010", file = f)
  expect_true(file.exists(f))
  expect_s3_class(x, "multi_longi")
  expect_equal(x$params$cluster_by, "joint")
  expect_equal(nlevels(x$cluster), 4)
  expect_s4_class(x$heatmap, "HeatmapList")

  expect_s3_class(multi_lineplot(x), "ggplot")
  expect_error(multi_cluster(multi_layers(), cluster_by = "joint"),
               "time column")
})

test_that("dual_heatmap matches the two-layer multi-omics result", {
  data("heatmap_example", package = "ModuleViz")
  layers <- multi_layers()
  obj <- longitudinal_cluster(layers$RNA, heatmap_example$meta,
                              sample_col = "SampleID", time_col = "TimeHr",
                              k = 4, stability = FALSE)

  d <- dual_heatmap(obj, layers$ATAC, file = tempfile(fileext = ".pdf"))
  m <- multi_cluster(list(RNA = obj, ATAC = layers$ATAC))
  expect_identical(d$z_secondary, m$layers$ATAC$z)
  expect_identical(d$genes, m$genes)
})
