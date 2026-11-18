test_that("the simulated single-cell data has the documented shape", {
  sce <- simulate_single_cell_example(cells_per_sample = 5)
  expect_s4_class(sce, "SummarizedExperiment")
  expect_equal(dim(sce), c(60L, 16L * 2L * 5L))
  expect_true(all(c("SampleID", "Condition", "TimeHr", "Replicate",
                    "CellType") %in%
                    names(SummarizedExperiment::colData(sce))))
  again <- simulate_single_cell_example(cells_per_sample = 5)
  expect_identical(SummarizedExperiment::assay(sce),
                   SummarizedExperiment::assay(again))
})

test_that("pseudobulk sums match a manual aggregation", {
  sce <- simulate_single_cell_example(cells_per_sample = 8)
  pb <- pseudobulk_cells(sce, time_col = "TimeHr", sample_col = "SampleID",
                         group_col = "Condition", min_cells = 1)
  expect_s4_class(pb, "SummarizedExperiment")
  expect_equal(ncol(pb), 16)
  expect_setequal(colnames(pb), unique(sce$SampleID))
  counts <- SummarizedExperiment::assay(sce, "counts")
  manual <- vapply(colnames(pb), function(s)
    rowSums(counts[, sce$SampleID == s, drop = FALSE]), numeric(nrow(counts)))
  expect_equal(SummarizedExperiment::assay(pb, "counts"), manual,
               ignore_attr = TRUE)
  expect_true(all(pb$n_cells == 16))
  expect_true(all(c("counts", "logcpm") %in%
                    SummarizedExperiment::assayNames(pb)))
  cpm <- sweep(manual, 2, colSums(manual), "/") * 1e6
  expect_equal(SummarizedExperiment::assay(pb, "logcpm"), log2(cpm + 1),
               ignore_attr = TRUE)
})

test_that("mean aggregation, sparse input, and matrix + meta agree", {
  sce <- simulate_single_cell_example(cells_per_sample = 6)
  counts <- SummarizedExperiment::assay(sce, "counts")
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  dense <- pseudobulk_cells(counts, meta = cd, time_col = "TimeHr",
                            sample_col = "SampleID", min_cells = 1)
  sparse <- pseudobulk_cells(methods::as(counts, "CsparseMatrix"), meta = cd,
                             time_col = "TimeHr", sample_col = "SampleID",
                             min_cells = 1)
  expect_equal(SummarizedExperiment::assay(dense, "counts"),
               SummarizedExperiment::assay(sparse, "counts"))
  avg <- pseudobulk_cells(counts, meta = cd, time_col = "TimeHr",
                          sample_col = "SampleID", fun = "mean",
                          min_cells = 1)
  expect_equal(SummarizedExperiment::assay(avg, 1),
               SummarizedExperiment::assay(dense, "counts") / 12)
  expect_error(pseudobulk_cells(counts, time_col = "TimeHr"), "meta")
  expect_error(pseudobulk_cells(counts, meta = cd, time_col = "Nope"),
               "not found")
  expect_error(pseudobulk_cells(counts, meta = cd[-1, ], time_col = "TimeHr"),
               "rows")
})

test_that("cell types, small groups, and pseudotime bins are handled", {
  sce <- simulate_single_cell_example(cells_per_sample = 6)
  by_type <- pseudobulk_cells(sce, time_col = "TimeHr",
                              sample_col = "SampleID",
                              celltype_col = "CellType", min_cells = 1)
  expect_equal(ncol(by_type), 32)
  expect_true("CellType" %in% names(SummarizedExperiment::colData(by_type)))

  keep <- !(sce$SampleID == sce$SampleID[1] & seq_len(ncol(sce)) > 3)
  expect_message(
    fewer <- pseudobulk_cells(sce[, keep], time_col = "TimeHr",
                              sample_col = "SampleID", min_cells = 5),
    "Dropping 1 pseudobulk")
  expect_equal(ncol(fewer), 15)
  expect_error(pseudobulk_cells(sce, time_col = "TimeHr", min_cells = 1e6),
               "Fewer than 2")

  sce$pseudotime <- sce$TimeHr + seq_len(ncol(sce)) / 1e4
  bins <- pseudobulk_cells(sce, time_col = "pseudotime", n_bins = 5)
  expect_equal(ncol(bins), 5)
  expect_equal(sort(bins$pseudotime), 1:5)
  expect_false(is.unsorted(bins$pseudotime_mean[order(bins$pseudotime)]))
  expect_error(pseudobulk_cells(sce, time_col = "CellType", n_bins = 3),
               "numeric")
})

test_that("pseudobulk profiles recover the bulk temporal modules", {
  sce <- simulate_single_cell_example(cells_per_sample = 20)
  pb <- pseudobulk_cells(sce, time_col = "TimeHr", sample_col = "SampleID",
                         group_col = "Condition")
  obj <- longitudinal_cluster(pb, time_col = "TimeHr",
                              group_col = "Condition", k = 4,
                              aggregate_replicates = TRUE, stability = FALSE)
  expect_s3_class(obj, "longi")
  ## the three simulated programmes (genes 1-15, 16-30, 31-45) should each
  ## sit almost entirely in one module, and in three different modules
  genes <- rownames(example_data()$mat)
  module_of <- setNames(obj$membership$Cluster, obj$membership$ID)
  blocks <- list(genes[1:15], genes[16:30], genes[31:45])
  dominant <- vapply(blocks, function(g)
    names(which.max(table(module_of[g]))), character(1))
  purity <- vapply(seq_along(blocks), function(i)
    mean(module_of[blocks[[i]]] == dominant[i]), numeric(1))
  expect_true(all(purity >= 0.8))
  expect_equal(anyDuplicated(dominant), 0)
  ## and the pseudobulk profiles track the bulk data they were simulated from
  bulk <- example_obj(aggregate_replicates = TRUE)
  expect_gt(cor(as.vector(obj$z[genes[1:45], ]),
                as.vector(bulk$z[genes[1:45], ])), 0.9)
  expect_s3_class(pattern_lineplot(obj), "ggplot")
})

test_that("Seurat objects are read when SeuratObject is available", {
  skip_if_not_installed("SeuratObject")
  sce <- simulate_single_cell_example(cells_per_sample = 5)
  counts <- methods::as(SummarizedExperiment::assay(sce, "counts"),
                        "CsparseMatrix")
  so <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = as.data.frame(SummarizedExperiment::colData(sce)))
  from_seurat <- pseudobulk_cells(so, time_col = "TimeHr",
                                  sample_col = "SampleID", min_cells = 1)
  from_sce <- pseudobulk_cells(sce, time_col = "TimeHr",
                               sample_col = "SampleID", min_cells = 1)
  expect_equal(SummarizedExperiment::assay(from_seurat, "counts"),
               SummarizedExperiment::assay(from_sce, "counts"))
  expect_error(longitudinal_cluster(so), "pseudobulk_cells")
})

test_that("variable-feature selection keeps the most variable rows", {
  ex <- example_data()
  top <- select_variable_features(ex$mat, n = 10)
  expect_equal(dim(top), c(10L, ncol(ex$mat)))
  v <- apply(ex$mat, 1, var)
  expect_setequal(rownames(top), names(sort(v, decreasing = TRUE))[1:10])
  expect_equal(nrow(select_variable_features(ex$mat, n = 1e6)), nrow(ex$mat))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = ex$mat))
  expect_equal(nrow(select_variable_features(se, n = 7)), 7)
  expect_error(select_variable_features(ex$mat, n = 0), "positive")
})
