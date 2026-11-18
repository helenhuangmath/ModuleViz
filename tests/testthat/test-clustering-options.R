test_that("every clustering method returns ordered, complete modules", {
  ex <- example_data()
  for (m in c("kmeans", "hierarchical", "pam", "clara", "fuzzy")) {
    obj <- longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr",
                                group_col = "Condition", k = 4, method = m,
                                stability = FALSE)
    expect_s3_class(obj, "longi")
    expect_equal(nlevels(obj$cluster), 4)
    expect_equal(sum(table(obj$cluster)), nrow(ex$mat))
    expect_setequal(rownames(obj$z), rownames(ex$mat))
    expect_equal(obj$params$method, m)
  }
  expect_error(longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr",
                                    method = "nope"))
})

test_that("fuzzy clustering reports a membership score", {
  ex <- example_data()
  obj <- longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr", k = 4,
                              method = "fuzzy", stability = FALSE)
  expect_named(obj$membership, c("ID", "Cluster", "MaxAbsZ", "Membership"))
  expect_true(all(obj$membership$Membership > 0 &
                    obj$membership$Membership <= 1))
  expect_error(longitudinal_cluster(ex$mat[1:8, ], ex$meta,
                                    time_col = "TimeHr", k = 6,
                                    method = "fuzzy", stability = FALSE),
               "fuzzy")
})

test_that("correlation distance tolerates constant features", {
  ex <- example_data()
  z <- zscore_rows(ex$mat)
  z[1, ] <- 0
  obj <- longitudinal_cluster(z, ex$meta, time_col = "TimeHr", k = 3,
                              method = "hierarchical",
                              dist_method = "correlation", is_zscore = TRUE,
                              stability = FALSE)
  expect_equal(nlevels(obj$cluster), 3)
})

test_that("clustering is reproducible and leaves the RNG stream alone", {
  ex <- example_data()
  set.seed(42)
  before <- runif(1)
  set.seed(42)
  a <- longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr", k = 4,
                            n_resample = 3)
  after <- runif(1)
  expect_equal(before, after)
  b <- longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr", k = 4,
                            n_resample = 3)
  expect_identical(a$membership, b$membership)
  expect_identical(a$stability, b$stability)
})

test_that("modules are ordered from early-peaking to late-peaking", {
  obj <- example_obj(aggregate_replicates = TRUE)
  cen <- as.matrix(obj$centroids[, -1])
  ## peak position is evaluated on the time-point means pooled over groups
  tp <- vapply(c("0", "6", "24", "48"), function(t)
    rowMeans(cen[, sub("^.*\\|", "", colnames(cen)) == t, drop = FALSE]),
    numeric(nrow(cen)))
  expect_false(is.unsorted(apply(tp, 1, which.max)))
})

test_that("time may be a factor or character and keeps its order", {
  ex <- example_data()
  meta <- ex$meta
  stage <- c("0" = "Naive", "6" = "Early", "24" = "Mid", "48" = "Late")
  meta$Stage <- factor(stage[as.character(meta$TimeHr)],
                       levels = c("Naive", "Early", "Mid", "Late"))
  num <- longitudinal_cluster(ex$mat, meta, time_col = "TimeHr",
                              group_col = "Condition", k = 4,
                              aggregate_replicates = TRUE, stability = FALSE)
  fac <- longitudinal_cluster(ex$mat, meta, time_col = "Stage",
                              group_col = "Condition", k = 4,
                              aggregate_replicates = TRUE, stability = FALSE)
  expect_equal(fac$params$time_levels, c("Naive", "Early", "Mid", "Late"))
  expect_equal(unname(fac$z), unname(num$z))
  expect_identical(fac$membership, num$membership)
  expect_false(anyNA(fac$meta$Stage))

  p <- pattern_lineplot(fac)
  expect_equal(levels(p$data$Time), c("Naive", "Early", "Mid", "Late"))
  expect_false(anyNA(p$data$Time))

  ## character time: order of first appearance
  meta$StageChr <- as.character(meta$Stage)
  meta <- meta[order(meta$TimeHr), ]
  chr <- longitudinal_cluster(ex$mat, meta, time_col = "StageChr", k = 4,
                              stability = FALSE)
  expect_equal(chr$params$time_levels, c("Naive", "Early", "Mid", "Late"))
})

test_that("aggregating replicates keeps constant metadata columns", {
  ex <- example_data()
  meta <- ex$meta
  meta$Batch <- "B1"
  meta$Arm <- paste0(meta$Condition, "_arm")
  obj <- longitudinal_cluster(ex$mat, meta, time_col = "TimeHr",
                              group_col = "Condition", k = 4,
                              aggregate_replicates = TRUE, stability = FALSE)
  expect_equal(ncol(obj$z), 8)
  expect_true(all(c("Batch", "Arm", "Condition", "TimeHr") %in%
                    names(obj$meta)))
  expect_false("Replicate" %in% names(obj$meta))
  expect_type(obj$meta$TimeHr, "double")
  expect_identical(colnames(obj$z), obj$meta$SampleID)
})

test_that("a SummarizedExperiment is accepted directly", {
  ex <- example_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = 2^ex$mat, logcounts = ex$mat),
    colData = ex$meta)
  colnames(se) <- ex$meta$SampleID
  from_se <- longitudinal_cluster(se, time_col = "TimeHr",
                                  group_col = "Condition", k = 4,
                                  stability = FALSE)
  from_mat <- example_obj()
  expect_identical(from_se$membership, from_mat$membership)
  expect_error(longitudinal_cluster(se, time_col = "TimeHr", assay = "nope"),
               "not found")
  expect_error(longitudinal_cluster(ex$mat, time_col = "TimeHr"), "meta")
})

test_that("input problems give clear errors", {
  ex <- example_data()
  expect_error(longitudinal_cluster(ex$mat, ex$meta, time_col = "Nope"),
               "not found")
  bad <- ex$meta
  bad$TimeHr[1] <- NA
  expect_error(longitudinal_cluster(ex$mat, bad, time_col = "TimeHr"),
               "missing")
  dup <- ex$meta
  dup$SampleID[2] <- dup$SampleID[1]
  expect_error(longitudinal_cluster(ex$mat, dup, time_col = "TimeHr"),
               "duplicated")
  with_na <- ex$mat
  with_na[1, 1] <- NA
  expect_message(
    obj <- longitudinal_cluster(with_na, ex$meta, time_col = "TimeHr", k = 3,
                                stability = FALSE),
    "Dropping 1 feature")
  expect_equal(nrow(obj$z), nrow(ex$mat) - 1)
  ## sample IDs that look like numbers must be matched by name
  m <- ex$mat
  colnames(m) <- as.character(100 + seq_len(ncol(m)))
  meta <- ex$meta
  meta$SampleID <- 100 + seq_len(ncol(m))
  obj <- longitudinal_cluster(m, meta, time_col = "TimeHr", k = 3,
                              stability = FALSE)
  expect_equal(ncol(obj$z), ncol(m))
})

test_that("choose_k copes with short ranges and other methods", {
  z <- zscore_rows(example_data()$mat)
  one <- choose_k(z, k_range = 3)
  expect_equal(attr(one, "suggested_k"), 3)
  two <- choose_k(z, k_range = 2:3, method = "pam")
  expect_true(attr(two, "suggested_k") %in% 2:3)
  expect_error(choose_k(z, k_range = 1000), "No usable k")
  f <- tempfile(fileext = ".pdf")
  choose_k(z, k_range = 2:5, file = f, chosen_k = 4)
  expect_true(file.exists(f))
})

test_that("labels match whole tokens before substrings", {
  ex <- example_data()
  mat <- ex$mat[1:6, ]
  rownames(mat) <- c("peak1_Il2", "peak2_Il2ra", "peak3_Il21", "peak4_Cd8a",
                     "peak5_Tox", "peak6_Tcf7")
  obj <- longitudinal_cluster(mat, ex$meta, time_col = "TimeHr", k = 2,
                              stability = FALSE)
  expect_equal(resolve_labels(obj, features = "Il2")$label, "peak1_Il2")
  expect_message(none <- resolve_labels(obj, features = "Nope"), "not found")
  expect_equal(nrow(none), 0)
})

test_that("trend plots default to one column and can be wrapped", {
  obj <- example_obj()
  one <- pattern_lineplot(obj)
  two <- pattern_lineplot(obj, ncol = 2)
  layout_of <- function(p) ggplot2::ggplot_build(p)$layout$layout
  expect_equal(max(layout_of(one)$COL), 1)
  expect_equal(max(layout_of(one)$ROW), 4)
  expect_equal(max(layout_of(two)$COL), 2)
  expect_equal(max(layout_of(two)$ROW), 2)

  ex <- example_data()
  ungrouped <- longitudinal_cluster(ex$mat, ex$meta, time_col = "TimeHr",
                                    k = 4, stability = FALSE)
  expect_equal(max(layout_of(pattern_lineplot(ungrouped))$COL), 1)
  f <- tempfile(fileext = ".pdf")
  pattern_lineplot(obj, ncol = 2, file = f)
  expect_true(file.exists(f))
})

test_that("heatmaps and correlation plots draw for every option", {
  obj <- example_obj(aggregate_replicates = TRUE)
  with_null_device({
    ht <- longitudinal_heatmap(obj, top_n_label = 1,
                               zscore_palette = "viridis",
                               show_column_names = FALSE)
    expect_s4_class(ht, "Heatmap")
    expect_warning(longitudinal_heatmap(obj, annotation_cols = "Nope"),
                   "not found")
    res <- module_correlation_heatmap(obj, clusters = "Cluster1",
                                      max_genes = 10, select = "random")
    expect_named(res, "Cluster1")
  })
  for (ext in c("png", "pdf", "svg")) {
    f <- tempfile(fileext = paste0(".", ext))
    longitudinal_heatmap(obj, file = f, width = 6, height = 6)
    expect_gt(file.size(f), 0)
  }
  summ <- module_correlation_summary(obj)
  expect_equal(nrow(summ), 4)
  prefix <- file.path(tempdir(), "mv_test")
  write_memberships(obj, prefix = prefix)
  expect_true(file.exists(paste0(prefix, "_memberships.txt")))
  expect_true(file.exists(paste0(prefix, "_centroids.txt")))
})
