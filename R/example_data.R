#' List supported public time-course example datasets.
#'
#' These examples come from the Bioconductor `timecoursedata` experiment-data
#' package. The returned keys can be passed to [load_real_timecourse_example()].
#'
#' @return data.table with dataset keys and short descriptions.
available_real_timecourse_datasets <- function() {
  data.table::data.table(
    dataset = c("sorghum_leaf", "sorghum_root", "influenza_mouse"),
    package_object = c("varoquaux2019leaf", "varoquaux2019root", "shoemaker2015"),
    assay = c("RNA-seq", "RNA-seq", "microarray"),
    description = c(
      "Sorghum leaf samples across weekly drought/control time courses.",
      "Sorghum root samples across weekly drought/control time courses.",
      "Mouse lung samples across influenza infection time points."
    )
  )
}

#' Load a real public time-course example in ModuleViz format.
#'
#' This helper downloads/loads one dataset from the Bioconductor
#' `timecoursedata` package and returns a small matrix plus matching metadata
#' ready for [longitudinal_cluster()]. The matrix is trimmed to the most variable
#' features by default so the examples run quickly on a laptop.
#'
#' @param dataset one of `"sorghum_leaf"`, `"sorghum_root"`,
#'   `"influenza_mouse"`.
#' @param top_n_features number of most-variable features to keep. Set `Inf` to
#'   keep all features.
#' @param genotype for sorghum datasets, keep one genotype. Use `NULL` to keep
#'   all genotypes.
#' @param condition for sorghum datasets, optional condition filter.
#' @param group for influenza, optional group filter.
#' @param min_week for sorghum, minimum week retained. The default removes week 2
#'   because it is control-only in the common tutorial subset.
#' @param log_transform if TRUE, use log2(x + 1) before feature selection.
#' @return list with `mat`, `meta`, `sample_col`, `time_col`, `group_col`,
#'   `source`, and `citation_hint`.
load_real_timecourse_example <- function(dataset = c("sorghum_leaf", "sorghum_root",
                                                     "influenza_mouse"),
                                         top_n_features = 1000,
                                         genotype = "BT642",
                                         condition = NULL,
                                         group = NULL,
                                         min_week = 3,
                                         log_transform = TRUE) {
  dataset <- match.arg(dataset)
  if (!requireNamespace("timecoursedata", quietly = TRUE)) {
    stop(
      "The real-data examples need the Bioconductor package 'timecoursedata'.\n",
      "Install it with:\n",
      "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('timecoursedata')",
      call. = FALSE
    )
  }

  object_name <- switch(dataset,
    sorghum_leaf = "varoquaux2019leaf",
    sorghum_root = "varoquaux2019root",
    influenza_mouse = "shoemaker2015"
  )
  env <- new.env(parent = emptyenv())
  utils::data(list = object_name, package = "timecoursedata", envir = env)
  raw <- get(object_name, envir = env)
  dat <- .extract_timecoursedata_list(raw)
  mat <- as.matrix(dat$data)
  mode(mat) <- "numeric"
  meta <- data.table::as.data.table(dat$meta)

  sample_ids <- .infer_sample_ids(mat, meta)
  meta[, SampleID := sample_ids]

  if (dataset %in% c("sorghum_leaf", "sorghum_root")) {
    if (!"Week" %in% names(meta) && "Timepoint" %in% names(meta)) {
      meta[, Week := as.numeric(Timepoint)]
    }
    if (!"Week" %in% names(meta)) {
      stop("Could not find a Week or Timepoint column in the sorghum metadata.")
    }
    if (!is.null(genotype) && "Genotype" %in% names(meta)) {
      meta <- meta[Genotype %in% genotype]
    }
    if (!is.null(condition) && "Condition" %in% names(meta)) {
      meta <- meta[Condition %in% condition]
    }
    if (!is.null(min_week)) {
      meta <- meta[as.numeric(Week) >= min_week]
    }
    time_col <- "Week"
    group_col <- if ("Condition" %in% names(meta)) "Condition" else NULL
    source <- "Bioconductor timecoursedata: Varoquaux 2019 sorghum drought time-course"
  } else {
    if (!"Timepoint" %in% names(meta)) {
      stop("Could not find a Timepoint column in the influenza metadata.")
    }
    if (!is.null(group) && "Group" %in% names(meta)) {
      meta <- meta[Group %in% group]
    }
    time_col <- "Timepoint"
    group_col <- if ("Group" %in% names(meta)) "Group" else NULL
    source <- "Bioconductor timecoursedata: Shoemaker 2015 influenza mouse lung time-course"
  }

  keep_samples <- intersect(meta$SampleID, colnames(mat))
  if (length(keep_samples) < 2) {
    stop("Fewer than two matching samples after filtering; relax the filters.")
  }
  meta <- meta[match(keep_samples, SampleID)]
  mat <- mat[, keep_samples, drop = FALSE]

  if (log_transform) {
    mat <- log2(mat + 1)
  }
  mat <- .keep_most_variable_features(mat, top_n_features)

  list(
    mat = mat,
    meta = as.data.frame(meta),
    sample_col = "SampleID",
    time_col = time_col,
    group_col = group_col,
    source = source,
    citation_hint = "Run citation('timecoursedata') and cite the original dataset publication."
  )
}

.extract_timecoursedata_list <- function(x) {
  if (is.list(x) && all(c("data", "meta") %in% names(x))) {
    return(list(data = x$data, meta = x$meta))
  }
  data_attr <- attr(x, "data")
  meta_attr <- attr(x, "meta")
  if (!is.null(data_attr) && !is.null(meta_attr)) {
    return(list(data = data_attr, meta = meta_attr))
  }
  stop("Unsupported timecoursedata object format.")
}

.infer_sample_ids <- function(mat, meta) {
  candidates <- c("SampleID", "sample_id", "Sample", "Barcode", "libraryName")
  for (col in intersect(candidates, names(meta))) {
    ids <- as.character(meta[[col]])
    if (sum(ids %in% colnames(mat)) >= 2) {
      return(ids)
    }
  }
  rn <- rownames(meta)
  if (!is.null(rn) && sum(rn %in% colnames(mat)) >= 2) {
    return(rn)
  }
  if (nrow(meta) == ncol(mat)) {
    return(colnames(mat))
  }
  stop("Could not infer sample IDs linking metadata rows to matrix columns.")
}

.keep_most_variable_features <- function(mat, top_n_features) {
  if (is.infinite(top_n_features) || top_n_features >= nrow(mat)) {
    return(mat)
  }
  vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(max(2L, top_n_features))]
  mat[keep, , drop = FALSE]
}
