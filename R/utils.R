## =============================================================================
## utils.R
##
## Small internal helpers shared by the clustering, plotting, and single-cell
## code: reproducible random-number use, input coercion from Bioconductor
## containers, and time-axis handling.
## =============================================================================

#' Evaluate an expression with a fixed RNG seed, then restore the RNG state.
#'
#' `seed = NULL` leaves the random-number generator untouched, so the caller's
#' own `set.seed()` (if any) governs the result.
#' @noRd
.moduleviz_with_seed <- function(seed, expr) {
    if (is.null(seed)) {
        return(expr)
    }
    withr::with_seed(seed, expr)
}

#' Ordered levels of a time column.
#'
#' Numeric time is sorted numerically, a factor keeps its level order, and
#' character time keeps the order of first appearance.
#' @noRd
.moduleviz_time_levels <- function(x) {
    if (is.factor(x)) {
        return(intersect(levels(x), as.character(x)))
    }
    if (is.numeric(x)) {
        return(as.character(sort(unique(x))))
    }
    unique(as.character(x[!is.na(x)]))
}

#' Names of the metadata columns that are constant within every group.
#'
#' Used when samples (or cells) are aggregated: a column that takes a single
#' value inside each aggregation unit can be carried over to the aggregated
#' metadata; anything else cannot.
#' @noRd
.moduleviz_constant_cols <- function(meta, key) {
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    is_const <- vapply(names(meta), function(cn) {
        v <- meta[[cn]]
        if (is.list(v)) {
            return(FALSE)
        }
        pairs <- unique(data.frame(k = key, v = as.character(v),
                                   stringsAsFactors = FALSE))
        !anyDuplicated(pairs$k)
    }, logical(1))
    names(meta)[is_const]
}

#' Pick the assay to analyse from a SummarizedExperiment.
#'
#' An explicit `assay` wins.  Otherwise the first log-scale assay found is
#' used, falling back to the first assay.
#' @noRd
.moduleviz_pick_assay <- function(x, assay = NULL,
                                  prefer = c("logcpm", "logcounts",
                                             "normcounts", "vst", "rlog")) {
    available <- SummarizedExperiment::assayNames(x)
    if (!is.null(assay)) {
        if (is.character(assay) && !assay %in% available) {
            stop("assay '", assay, "' not found. Available: ",
                 paste(available, collapse = ", "), call. = FALSE)
        }
        return(assay)
    }
    hit <- intersect(prefer, available)
    if (length(hit)) {
        return(hit[1])
    }
    if (length(available) && identical(available[1], "counts")) {
        message("Using the raw 'counts' assay; ModuleViz expects ",
                "log-normalised values. Set `assay` to choose another.")
    }
    1L
}

#' Turn any supported input into a features x samples matrix plus metadata.
#'
#' Accepts a plain matrix / data.frame (returned untouched together with
#' `meta`) or a SummarizedExperiment (including SingleCellExperiment), whose
#' `colData` supplies the sample metadata when `meta` is not given.
#' @return list(mat, meta, sample_col)
#' @noRd
.moduleviz_unpack <- function(x, meta = NULL, sample_col = "SampleID",
                              assay = NULL) {
    if (inherits(x, "Seurat")) {
        stop("`mat` is a Seurat object. Summarise the cells first with ",
             "pseudobulk_cells(), then pass the result here.", call. = FALSE)
    }
    if (!methods::is(x, "SummarizedExperiment")) {
        return(list(mat = x, meta = meta, sample_col = sample_col))
    }
    if (is.null(colnames(x))) {
        stop("The SummarizedExperiment has no column names.", call. = FALSE)
    }
    assay <- .moduleviz_pick_assay(x, assay)
    m <- as.matrix(SummarizedExperiment::assay(x, assay))
    if (is.null(rownames(m))) {
        stop("The SummarizedExperiment has no row names.", call. = FALSE)
    }
    if (is.null(meta)) {
        meta <- as.data.frame(SummarizedExperiment::colData(x),
                              stringsAsFactors = FALSE)
        ## the column names are the sample IDs; never overwrite a metadata
        ## column that happens to share the name but holds something else
        if (sample_col %in% names(meta) &&
            !identical(as.character(meta[[sample_col]]), colnames(x))) {
            sample_col <- "ColumnID"
        }
        meta[[sample_col]] <- colnames(x)
        rownames(meta) <- NULL
    }
    list(mat = m, meta = meta, sample_col = sample_col)
}
