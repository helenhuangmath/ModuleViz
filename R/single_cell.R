## =============================================================================
## single_cell.R
##
## Entry points for single-cell data.  ModuleViz clusters features across an
## ordered set of samples, so cells are first summarised into pseudobulk
## profiles (one per sample / cell type / condition / time point, or one per
## pseudotime bin); the result then runs through the same clustering and
## plotting functions as bulk data.
## =============================================================================

#' Pull a features x cells matrix and per-cell metadata out of a container.
#' @noRd
.moduleviz_cells <- function(x, meta = NULL, assay = NULL, fun = "sum") {
    if (inherits(x, "Seurat")) {
        if (!requireNamespace("SeuratObject", quietly = TRUE)) {
            stop("Reading a Seurat object needs the 'SeuratObject' package.",
                 call. = FALSE)
        }
        layer <- if (fun == "sum") "counts" else "data"
        if (is.null(assay)) {
            assay <- SeuratObject::DefaultAssay(x)
        }
        ns <- asNamespace("SeuratObject")
        mat <- if (exists("LayerData", envir = ns, inherits = FALSE)) {
            tryCatch(
                SeuratObject::LayerData(x, assay = assay, layer = layer),
                error = function(e) {
                    stop("Could not read layer '", layer, "' of assay '",
                         assay, "'. If the assay is split into several ",
                         "layers, run JoinLayers() first. (",
                         conditionMessage(e), ")", call. = FALSE)
                })
        } else {
            SeuratObject::GetAssayData(x, assay = assay, slot = layer)
        }
        if (is.null(meta)) {
            meta <- methods::slot(x, "meta.data")
        }
        source_assay <- if (fun == "sum") "counts" else "logcounts"
    } else if (methods::is(x, "SummarizedExperiment")) {
        available <- SummarizedExperiment::assayNames(x)
        if (is.null(assay)) {
            want <- if (fun == "sum") "counts" else "logcounts"
            assay <- if (want %in% available) want else 1L
        } else if (is.character(assay) && !assay %in% available) {
            stop("assay '", assay, "' not found. Available: ",
                 paste(available, collapse = ", "), call. = FALSE)
        }
        mat <- SummarizedExperiment::assay(x, assay)
        if (is.null(meta)) {
            meta <- as.data.frame(SummarizedExperiment::colData(x),
                                  stringsAsFactors = FALSE)
        }
        source_assay <- if (is.character(assay)) assay else available[assay]
        if (is.null(source_assay) || is.na(source_assay)) {
            source_assay <- if (fun == "sum") "counts" else "mean"
        }
    } else {
        mat <- x
        if (is.null(meta)) {
            stop("`meta` (one row per cell) is required when `x` is a matrix.",
                 call. = FALSE)
        }
        source_assay <- if (fun == "sum") "counts" else "mean"
    }

    if (!(is.matrix(mat) || methods::is(mat, "Matrix"))) {
        mat <- as.matrix(mat)
    }
    if (is.null(rownames(mat))) {
        stop("The expression matrix has no feature (row) names.",
             call. = FALSE)
    }
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (nrow(meta) != ncol(mat)) {
        stop("`meta` has ", nrow(meta), " rows but the matrix has ",
             ncol(mat), " cells.", call. = FALSE)
    }
    ## align by cell name when both sides carry one
    cell_names <- colnames(mat)
    if (!is.null(cell_names) && !is.null(rownames(meta)) &&
        !identical(rownames(meta), cell_names) &&
        all(cell_names %in% rownames(meta))) {
        meta <- meta[cell_names, , drop = FALSE]
    }
    list(mat = mat, meta = meta, source_assay = source_assay)
}

#' Summarise single cells into pseudobulk profiles.
#'
#' Cells are grouped by every metadata column supplied (`sample_col`,
#' `celltype_col`, `group_col`, `time_col`) and their expression is summed or
#' averaged within each group.  The result is a `SummarizedExperiment` with one
#' column per pseudobulk profile, ready for [longitudinal_cluster()].
#'
#' With `fun = "sum"` (the default, for raw counts) the summed counts are also
#' converted to log2 counts-per-million, stored in the `logcpm` assay.  With
#' `fun = "mean"` the values are averaged as they are, which suits data that
#' are already log-normalised.
#'
#' For data without discrete time points, give a continuous ordering such as
#' pseudotime as `time_col` and set `n_bins`: cells are split into `n_bins`
#' equal-sized bins along that ordering and the bin index becomes the time.
#'
#' @param x a `SingleCellExperiment` / `SummarizedExperiment`, a `Seurat`
#'   object, or a features x cells matrix (dense or sparse).
#' @param meta per-cell metadata, one row per cell.  Required for a matrix;
#'   optional otherwise, where it replaces the object's own cell metadata.
#' @param time_col metadata column giving the time point (or pseudotime) of
#'   each cell.
#' @param sample_col optional column identifying the biological sample /
#'   replicate each cell came from.
#' @param group_col optional column for a condition, genotype, or treatment.
#' @param celltype_col optional column of cell-type labels; when given, every
#'   cell type gets its own pseudobulk profiles.
#' @param assay assay (or Seurat assay) to read.  Defaults to `counts` for
#'   `fun = "sum"` and `logcounts` for `fun = "mean"`.
#' @param fun `"sum"` or `"mean"`.
#' @param normalize add a log2-CPM assay.  Defaults to `TRUE` for
#'   `fun = "sum"` and is ignored for `fun = "mean"`.
#' @param min_cells pseudobulk profiles built from fewer cells are dropped.
#' @param n_bins optional number of equal-sized bins to cut a continuous
#'   `time_col` (for example pseudotime) into.
#' @param features optional subset of features to keep.
#' @param prior_count,scale_factor log2-CPM settings:
#'   `log2(count / library_size * scale_factor + prior_count)`.
#' @return A `SummarizedExperiment`.  Its `colData` holds the grouping
#'   columns, every other metadata column that is constant within a profile,
#'   and `n_cells`.
#' @examples
#' sce <- simulate_single_cell_example(cells_per_sample = 10)
#' pb <- pseudobulk_cells(sce, time_col = "TimeHr", sample_col = "SampleID",
#'                        group_col = "Condition", min_cells = 5)
#' pb
#' obj <- longitudinal_cluster(pb, time_col = "TimeHr",
#'                             group_col = "Condition", k = 4,
#'                             stability = FALSE)
#' obj
pseudobulk_cells <- function(x, meta = NULL, time_col,
                             sample_col = NULL, group_col = NULL,
                             celltype_col = NULL, assay = NULL,
                             fun = c("sum", "mean"), normalize = NULL,
                             min_cells = 10, n_bins = NULL, features = NULL,
                             prior_count = 1, scale_factor = 1e6) {
    fun <- match.arg(fun)
    if (missing(time_col) || is.null(time_col)) {
        stop("`time_col` is required.", call. = FALSE)
    }
    if (is.null(normalize)) {
        normalize <- fun == "sum"
    }
    cells <- .moduleviz_cells(x, meta = meta, assay = assay, fun = fun)
    mat <- cells$mat
    cd <- cells$meta

    by <- unique(c(sample_col, celltype_col, group_col, time_col))
    absent <- setdiff(by, names(cd))
    if (length(absent)) {
        stop("Column(s) not found in the cell metadata: ",
             paste(absent, collapse = ", "), call. = FALSE)
    }
    if (!is.null(features)) {
        features <- intersect(features, rownames(mat))
        if (!length(features)) {
            stop("None of `features` are in the expression matrix.",
                 call. = FALSE)
        }
        mat <- mat[features, , drop = FALSE]
    }

    ## ---- drop cells with missing grouping values ---------------------------
    ok <- stats::complete.cases(cd[, by, drop = FALSE])
    if (!all(ok)) {
        message("Dropping ", sum(!ok),
                " cell(s) with missing values in the grouping columns.")
        mat <- mat[, ok, drop = FALSE]
        cd <- cd[ok, , drop = FALSE]
    }

    ## ---- optional binning of a continuous ordering -------------------------
    if (!is.null(n_bins)) {
        tv <- cd[[time_col]]
        if (!is.numeric(tv)) {
            stop("`n_bins` needs a numeric `time_col` (e.g. pseudotime).",
                 call. = FALSE)
        }
        n_bins <- as.integer(n_bins)
        if (is.na(n_bins) || n_bins < 2) {
            stop("`n_bins` must be an integer >= 2.", call. = FALSE)
        }
        bin <- cut(rank(tv, ties.method = "first"), breaks = n_bins,
                   labels = FALSE)
        cd[[paste0(time_col, "_mean")]] <- stats::ave(tv, bin)
        cd[[time_col]] <- bin
    }

    ## ---- pseudobulk units --------------------------------------------------
    ## Name each unit from the columns that actually tell units apart: a
    ## column already implied by the ones before it (the condition of a
    ## sample, say) would only make the names longer.
    key <- NULL
    for (cn in by) {
        candidate <- if (is.null(key)) as.character(cd[[cn]])
                     else paste(key, cd[[cn]], sep = "_")
        if (is.null(key) || length(unique(candidate)) > length(unique(key))) {
            key <- candidate
        }
    }
    n_cells <- table(key)
    keep_units <- unique(key)[n_cells[unique(key)] >= min_cells]
    n_small <- length(n_cells) - length(keep_units)
    if (n_small) {
        message("Dropping ", n_small, " pseudobulk profile(s) with fewer ",
                "than ", min_cells, " cells.")
    }
    if (length(keep_units) < 2) {
        stop("Fewer than 2 pseudobulk profiles have at least ", min_cells,
             " cells; lower `min_cells` or group the cells more coarsely.",
             call. = FALSE)
    }
    in_unit <- key %in% keep_units
    mat <- mat[, in_unit, drop = FALSE]
    cd <- cd[in_unit, , drop = FALSE]
    key <- key[in_unit]
    unit <- factor(key, levels = keep_units)

    ## ---- aggregate: (features x cells) %*% (cells x units) -----------------
    indicator <- Matrix::sparseMatrix(i = seq_along(unit),
                                      j = as.integer(unit), x = 1,
                                      dims = c(length(unit), nlevels(unit)))
    agg <- as.matrix(mat %*% indicator)
    dimnames(agg) <- list(rownames(mat), keep_units)
    unit_n <- as.integer(n_cells[keep_units])
    if (fun == "mean") {
        agg <- sweep(agg, 2, unit_n, "/")
    }

    ## ---- metadata: everything that is constant within a unit ---------------
    const <- .moduleviz_constant_cols(cd, key)
    col_data <- cd[match(keep_units, key), const, drop = FALSE]
    col_data$n_cells <- unit_n
    rownames(col_data) <- keep_units

    assays <- if (fun == "sum") {
        out <- list(counts = agg)
        if (isTRUE(normalize)) {
            lib <- colSums(agg)
            lib[lib == 0] <- 1
            out$logcpm <- log2(sweep(agg, 2, lib, "/") * scale_factor +
                                   prior_count)
        }
        out
    } else {
        stats::setNames(list(agg), cells$source_assay)
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = assays, colData = S4Vectors::DataFrame(col_data,
                                                        check.names = FALSE))
}

#' Keep the most variable features.
#'
#' Clustering every gene of a single-cell or bulk dataset is rarely useful:
#' most features are flat across the time course.  This keeps the `n` features
#' with the largest variance across samples.
#'
#' @param x a features x samples matrix or a `SummarizedExperiment`.
#' @param n number of features to keep.
#' @param assay assay used to rank the features when `x` is a
#'   `SummarizedExperiment`; chosen as in [longitudinal_cluster()] by default.
#' @return `x` restricted to its `n` most variable features, most variable
#'   first.
#' @examples
#' data("heatmap_example", package = "ModuleViz")
#' top <- select_variable_features(heatmap_example$mat, n = 20)
#' dim(top)
select_variable_features <- function(x, n = 2000, assay = NULL) {
    if (methods::is(x, "SummarizedExperiment")) {
        m <- as.matrix(SummarizedExperiment::assay(
            x, .moduleviz_pick_assay(x, assay)))
    } else {
        m <- as.matrix(x)
    }
    if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
        stop("`n` must be a single positive number.", call. = FALSE)
    }
    centred <- m - rowMeans(m, na.rm = TRUE)
    v <- rowMeans(centred^2, na.rm = TRUE)
    v[is.na(v)] <- 0
    keep <- order(v, decreasing = TRUE)[seq_len(min(n, nrow(m)))]
    x[keep, , drop = FALSE]
}

#' Simulate a small single-cell time course.
#'
#' Generates single-cell counts whose sample-level averages follow the temporal
#' patterns of the bundled [heatmap_example] data: 60 genes, two conditions,
#' four time points, and two replicates, with several cell types per sample.
#' Intended for examples, tests, and tutorials.
#'
#' @param cells_per_sample number of cells per sample and cell type.
#' @param cell_types names of the simulated cell types.
#' @param seed RNG seed; `NULL` uses the current RNG state.
#' @return A `SingleCellExperiment` (a `SummarizedExperiment` if the
#'   SingleCellExperiment package is not installed) with a `counts` assay and
#'   the cell metadata columns `SampleID`, `Condition`, `TimeHr`, `Replicate`,
#'   and `CellType`.
#' @examples
#' sce <- simulate_single_cell_example(cells_per_sample = 10)
#' sce
simulate_single_cell_example <- function(cells_per_sample = 20,
                                         cell_types = c("T_cell", "Monocyte"),
                                         seed = 1) {
    env <- new.env(parent = emptyenv())
    utils::data("heatmap_example", package = "ModuleViz", envir = env)
    bulk <- env$heatmap_example$mat
    bulk_meta <- env$heatmap_example$meta
    cells_per_sample <- as.integer(cells_per_sample)
    if (is.na(cells_per_sample) || cells_per_sample < 1) {
        stop("`cells_per_sample` must be a positive integer.", call. = FALSE)
    }

    grid <- expand.grid(cell = seq_len(cells_per_sample),
                        CellType = cell_types,
                        sample = seq_len(nrow(bulk_meta)),
                        stringsAsFactors = FALSE)
    cell_meta <- data.frame(
        SampleID = bulk_meta$SampleID[grid$sample],
        Condition = bulk_meta$Condition[grid$sample],
        TimeHr = bulk_meta$TimeHr[grid$sample],
        Replicate = bulk_meta$Replicate[grid$sample],
        CellType = grid$CellType,
        stringsAsFactors = FALSE)
    cell_id <- sprintf("%s_%s_cell%03d", cell_meta$SampleID,
                       cell_meta$CellType, grid$cell)
    rownames(cell_meta) <- cell_id

    ## expected counts: the bulk pattern on a log2 scale around ~20 counts
    centred <- (bulk - mean(bulk)) / stats::sd(bulk)
    mu <- 20 * 2^centred[, grid$sample, drop = FALSE]
    type_shift <- stats::setNames(seq(0.8, 1.2, length.out = length(cell_types)),
                                  cell_types)

    counts <- .moduleviz_with_seed(seed, {
        lib <- stats::rlnorm(nrow(cell_meta), meanlog = 0, sdlog = 0.2)
        expected <- sweep(mu, 2, lib * type_shift[cell_meta$CellType], "*")
        matrix(stats::rnbinom(length(expected), mu = expected, size = 5),
               nrow = nrow(expected))
    })
    dimnames(counts) <- list(rownames(bulk), cell_id)

    col_data <- S4Vectors::DataFrame(cell_meta, check.names = FALSE)
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts), colData = col_data)
    } else {
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = counts), colData = col_data)
    }
}
