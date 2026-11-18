example_sets <- function() {
  g <- rownames(example_data()$mat)
  list(Early = g[1:15], Mid = g[16:30], Late = g[31:45],
       Mixed = g[seq(2, 60, by = 4)])
}

test_that("module enrichment reproduces the hypergeometric test", {
  obj <- example_obj()
  sets <- example_sets()
  enr <- module_enrichment(obj, sets, min_size = 3, min_overlap = 1)
  expect_s3_class(enr, "data.table")
  expect_named(enr, c("Cluster", "Term", "Overlap", "ModuleSize", "SetSize",
                      "Universe", "FoldEnrichment", "PValue", "Padj",
                      "Genes"))
  row <- enr[1]
  mod <- rownames(obj$z)[obj$cluster == row$Cluster]
  set <- sets[[row$Term]]
  shared <- intersect(mod, set)
  expect_equal(row$Overlap, length(shared))
  expect_equal(row$PValue,
               stats::phyper(length(shared) - 1, length(set), 60 - length(set),
                             length(mod), lower.tail = FALSE))
  expect_equal(enr$Padj, stats::p.adjust(enr$PValue, "BH"))
  expect_setequal(strsplit(row$Genes, ",")[[1]], shared)
  ## each simulated programme is the top hit of exactly one module
  top <- enr[, .SD[1], by = Cluster]
  expect_true(all(c("Early", "Mid", "Late") %in% top$Term))
})

test_that("gene sets can be a data frame or a GMT file", {
  obj <- example_obj()
  sets <- example_sets()
  from_list <- module_enrichment(obj, sets, min_size = 3)
  df <- data.frame(term = rep(names(sets), lengths(sets)),
                   gene = unlist(sets))
  expect_equal(module_enrichment(obj, df, min_size = 3), from_list)
  gmt <- tempfile(fileext = ".gmt")
  writeLines(vapply(names(sets), function(n)
    paste(c(n, "desc", sets[[n]]), collapse = "\t"), character(1)), gmt)
  expect_equal(module_enrichment(obj, gmt, min_size = 3), from_list)
  expect_error(module_enrichment(obj, list(1:3)), "named list")
  expect_error(module_enrichment(obj, "missing_file.gmt"), "not found")
})

test_that("universe, size filters, and gene mapping are respected", {
  obj <- example_obj()
  sets <- example_sets()
  wide <- module_enrichment(obj, sets, min_size = 3,
                            universe = c(rownames(obj$z),
                                         paste0("Other", 1:940)))
  expect_true(all(wide$Universe == 1000))
  narrow <- module_enrichment(obj, sets, min_size = 3)
  expect_lt(min(wide$PValue), min(narrow$PValue))
  expect_warning(none <- module_enrichment(obj, sets, min_size = 100),
                 "nothing was tested")
  expect_equal(nrow(none), 0)

  ex <- example_data()
  peaks <- ex$mat
  rownames(peaks) <- paste0("peak", seq_len(nrow(peaks)), "_",
                            rownames(peaks))
  pobj <- longitudinal_cluster(peaks, ex$meta, time_col = "TimeHr",
                               group_col = "Condition", k = 4,
                               stability = FALSE)
  mapped <- module_enrichment(pobj, sets, min_size = 3,
                              gene_of = function(id) sub("^.*_", "", id))
  expect_equal(mapped$PValue, narrow$PValue)
})

test_that("enrichment is drawn beside heatmaps and as a dot plot", {
  obj <- example_obj()
  enr <- module_enrichment(obj, example_sets(), min_size = 3)
  labels <- ModuleViz:::.moduleviz_enrichment_labels(enr, levels(obj$cluster),
                                                    top_n = 1)
  expect_length(labels, 4)
  expect_setequal(labels[nzchar(labels)], c("Early", "Mid", "Late"))

  f <- tempfile(fileext = ".pdf")
  ht <- longitudinal_heatmap(obj, enrichment = enr, top_n_label = 1,
                             file = f)
  expect_true(file.exists(f))
  expect_s4_class(ht, "Heatmap")
  ## nothing significant: the annotation is dropped, the heatmap still draws
  expect_message(
    longitudinal_heatmap(obj, enrichment = enr, enrichment_padj = 1e-300,
                         file = tempfile(fileext = ".pdf")),
    "omitted")
  expect_error(longitudinal_heatmap(obj, enrichment = data.frame(x = 1),
                                    file = tempfile(fileext = ".pdf")),
               "Cluster and Term")

  p <- enrichment_dotplot(enr)
  expect_s3_class(p, "ggplot")
  expect_error(enrichment_dotplot(enr, padj_cutoff = 1e-300), "No term")

  ex <- example_data()
  atac <- ex$mat[rep(rownames(ex$mat), each = 2), ]
  rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
  x <- multi_cluster(list(RNA = ex$mat, ATAC = atac), meta = ex$meta,
                     time_col = "TimeHr", group_col = "Condition", k = 4,
                     stability = FALSE)
  menr <- module_enrichment(x, example_sets(), min_size = 3)
  expect_equal(menr$PValue, enr$PValue)
  f2 <- tempfile(fileext = ".pdf")
  multi_heatmap(x, enrichment = menr, file = f2)
  expect_true(file.exists(f2))
})

test_that("multi-layer trend plots default to one column", {
  ex <- example_data()
  atac <- ex$mat[rep(rownames(ex$mat), each = 2), ]
  rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
  x <- multi_cluster(list(RNA = ex$mat, ATAC = atac), meta = ex$meta,
                     time_col = "TimeHr", group_col = "Condition", k = 4,
                     stability = FALSE)
  layout_of <- function(p) ggplot2::ggplot_build(p)$layout$layout
  expect_equal(max(layout_of(multi_lineplot(x))$COL), 1)
  expect_equal(max(layout_of(multi_lineplot(x, ncol = 2))$COL), 2)
  split <- layout_of(multi_lineplot(x, split_groups = TRUE))
  expect_equal(max(split$COL), 2)
  expect_equal(max(split$ROW), 4)
})
