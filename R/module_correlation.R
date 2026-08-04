## =============================================================================
## module_correlation.R
##
## Gene-gene correlation heatmaps, one per temporal module.
##
## longitudinal_cluster() tells you WHICH genes share a trajectory; it does not
## tell you how tightly they agree.  A module can look clean in the main heatmap
## and still be an average of two anti-correlated sub-programmes.  Correlating
## every gene in a module against every other gene, across the ordered samples,
## exposes that: a coherent module is a solid block of high correlation, while a
## module with sub-structure breaks into visible blocks.
## =============================================================================

#' Gene-gene correlation heatmap for each temporal module.
#'
#' For every module in a `longi` object, correlate its features against each
#' other across the (ordered) samples and draw the correlation matrix as a
#' heatmap.  Useful as a quality check on module coherence and for spotting
#' sub-programmes hiding inside a module.
#'
#' Large modules are subsampled to `max_genes` features, because a correlation
#' heatmap of several thousand genes is neither legible nor quick to render.
#' `select = "top"` keeps the strongest-signal features (largest max |z|), which
#' is usually what you want; `select = "random"` gives an unbiased view of the
#' whole module.
#'
#' @param obj a "longi" object from [longitudinal_cluster()].
#' @param clusters character vector of module names to draw (default: all).
#' @param method correlation method, "pearson" (default) or "spearman".
#' @param max_genes maximum features per module; larger modules are subsampled.
#' @param select how to subsample: "top" (strongest signal) or "random".
#' @param cluster_genes if TRUE (default), order genes within the heatmap by
#'   hierarchical clustering on correlation distance, revealing sub-blocks.
#' @param corr_palette continuous palette: "rdbu" (default) or "viridis".
#' @param corr_limits length-2 correlation range for the colour scale.
#' @param col_fun optional explicit colour function; overrides `corr_palette`.
#' @param show_gene_names show feature names; default TRUE when a module is
#'   small enough (<= 60 features shown) and FALSE otherwise.
#' @param label_fontsize font size for feature names.
#' @param legend_fontsize,title_fontsize font sizes for the legend and title.
#' @param file optional output path.  A `.pdf` path collects every module as one
#'   multi-page PDF; any other extension writes one file per module, with the
#'   module name inserted before the extension.
#' @param width,height figure size in inches.
#' @param title optional title prefix; the module name and its mean correlation
#'   are appended.
#' @param seed RNG seed used when `select = "random"`.
#' @return invisibly, a named list with one entry per module, each holding
#'   `cor` (the correlation matrix), `mean_cor` (mean off-diagonal correlation),
#'   `n_total`, `n_shown`, and `heatmap`.
module_correlation_heatmap <- function(obj,
                                       clusters = NULL,
                                       method = c("pearson", "spearman"),
                                       max_genes = Inf,
                                       select = c("top", "random"),
                                       cluster_genes = TRUE,
                                       corr_palette = "rdbu",
                                       corr_limits = c(-1, 1),
                                       col_fun = NULL,
                                       show_gene_names = NULL,
                                       label_top_n = 30,
                                       label_fontsize = 8,
                                       legend_fontsize = 10,
                                       title_fontsize = 11,
                                       file = NULL, width = 7, height = 6.5,
                                       title = NULL, seed = 1) {
  stopifnot(inherits(obj, "longi"))
  method <- match.arg(method)
  select <- match.arg(select)
  if (is.null(col_fun))
    col_fun <- .moduleviz_zscore_col_fun(corr_palette, corr_limits)

  z  <- obj$z
  cl <- obj$cluster
  if (ncol(z) < 3)
    stop("Need at least 3 samples to compute gene-gene correlations.")

  levs <- if (is.null(clusters)) levels(cl) else intersect(clusters, levels(cl))
  if (!length(levs))
    stop("None of the requested clusters exist in this object.")

  ## A .pdf destination holds every module as separate pages in one device;
  ## anything else gets one file per module.
  multipage <- !is.null(file) &&
    identical(tolower(sub("^.*\\.([^.]+)$", "\\1", file)), "pdf")
  if (multipage) {
    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  cluster_cols <- .moduleviz_cluster_cols(levels(cl), "paired")
  out <- list()

  for (lev in levs) {
    idx  <- which(cl == lev)
    subz <- z[idx, , drop = FALSE]
    n_total <- nrow(subz)

    ## constant rows have undefined correlation - drop them up front
    subz <- subz[apply(subz, 1, stats::sd) > 0, , drop = FALSE]
    if (nrow(subz) < 3) {
      message("skipping ", lev, ": fewer than 3 non-constant features")
      next
    }

    ## By default every feature is kept; max_genes is an optional cap.
    if (is.finite(max_genes) && nrow(subz) > max_genes) {
      if (select == "top") {
        keep <- order(apply(abs(subz), 1, max), decreasing = TRUE)[seq_len(max_genes)]
        subz <- subz[sort(keep), , drop = FALSE]
      } else {
        set.seed(seed)
        subz <- subz[sort(sample.int(nrow(subz), max_genes)), , drop = FALSE]
      }
    }

    cm <- stats::cor(t(subz), method = method)
    mean_cor <- mean(cm[lower.tri(cm)], na.rm = TRUE)
    n_shown  <- nrow(cm)

    ## Order rows/columns on the correlation matrix itself.  Clustering the
    ## correlation profiles instead would recompute an n x n correlation of an
    ## n x n matrix, which is needlessly slow once a module has 1000 features.
    ord <- if (cluster_genes) stats::hclust(stats::as.dist(1 - cm), "average") else FALSE

    ## Labelling.  Printing every row name is unreadable beyond a few dozen
    ## features, so past that we mark only the strongest `label_top_n` - the
    ## features with the largest |z|, i.e. those moving most across the series.
    show_names <- if (is.null(show_gene_names)) n_shown <= 60 else isTRUE(show_gene_names)
    right_anno <- NULL
    marked <- character(0)
    if (!show_names && label_top_n > 0) {
      strength <- apply(abs(subz), 1, max)
      at <- sort(order(strength, decreasing = TRUE)[seq_len(min(label_top_n, n_shown))])
      marked <- rownames(cm)[at]
      right_anno <- rowAnnotation(mark = anno_mark(
        at = at, labels = marked,
        labels_gp = gpar(fontsize = label_fontsize, fontface = "italic"),
        link_width = unit(6, "mm")))
    } else if (show_names) {
      marked <- rownames(cm)
    }

    ht_title <- sprintf("%s%s  (n=%d%s, mean r=%.2f)",
                        if (is.null(title)) "" else paste0(title, " - "),
                        lev, n_total,
                        if (n_shown < n_total) sprintf(", showing %d", n_shown) else "",
                        mean_cor)

    ht <- Heatmap(
      cm,
      name = sprintf("%s r", method),
      col = col_fun,
      column_title = ht_title,
      column_title_gp = gpar(fontsize = title_fontsize),
      cluster_rows = ord, cluster_columns = ord,
      show_row_names = show_names, show_column_names = show_names,
      row_names_gp = gpar(fontsize = label_fontsize),
      column_names_gp = gpar(fontsize = label_fontsize),
      show_row_dend = cluster_genes, show_column_dend = FALSE,
      left_annotation = rowAnnotation(
        Module = rep(lev, n_shown),
        col = list(Module = cluster_cols[lev]),
        show_annotation_name = FALSE, show_legend = FALSE,
        width = unit(3, "mm")),
      right_annotation = right_anno,
      use_raster = TRUE, raster_quality = 2,
      heatmap_legend_param = c(list(direction = "horizontal",
                                    title_position = "topcenter",
                                    legend_width = unit(3.5, "cm")),
                               .moduleviz_legend_gp(legend_fontsize)))

    draw_it <- function() draw(ht, heatmap_legend_side = "bottom")
    if (multipage) {
      draw_it()
    } else if (!is.null(file)) {
      per <- sub("(\\.[^.]+)$", paste0("_", lev, "\\1"), file)
      .moduleviz_write_plot(per, width, height, draw_it)
    } else {
      draw_it()
    }

    out[[lev]] <- list(cor = cm, mean_cor = mean_cor, n_total = n_total,
                       n_shown = n_shown, labelled = marked, heatmap = ht)
  }

  if (multipage) message("wrote ", file)
  invisible(out)
}

#' Summarise how internally coherent each module is.
#'
#' Convenience wrapper returning just the mean within-module correlation per
#' module, without drawing anything.  Low values flag modules that are averaging
#' over genes which do not actually agree.
#'
#' @inheritParams module_correlation_heatmap
#' @return data.table with Cluster, N, NShown, MeanCor.
module_correlation_summary <- function(obj, method = c("pearson", "spearman"),
                                       max_genes = Inf, select = c("top", "random"),
                                       seed = 1) {
  stopifnot(inherits(obj, "longi"))
  method <- match.arg(method); select <- match.arg(select)
  z <- obj$z; cl <- obj$cluster

  rbindlist(lapply(levels(cl), function(lev) {
    sub <- z[cl == lev, , drop = FALSE]
    n_total <- nrow(sub)
    sub <- sub[apply(sub, 1, stats::sd) > 0, , drop = FALSE]
    if (nrow(sub) < 3)
      return(data.table(Cluster = lev, N = n_total, NShown = nrow(sub),
                        MeanCor = NA_real_))
    if (is.finite(max_genes) && nrow(sub) > max_genes) {
      if (select == "top") {
        sub <- sub[order(apply(abs(sub), 1, max), decreasing = TRUE)[seq_len(max_genes)], ,
                   drop = FALSE]
      } else {
        set.seed(seed)
        sub <- sub[sample.int(nrow(sub), max_genes), , drop = FALSE]
      }
    }
    cm <- stats::cor(t(sub), method = method)
    data.table(Cluster = lev, N = n_total, NShown = nrow(sub),
               MeanCor = round(mean(cm[lower.tri(cm)], na.rm = TRUE), 3))
  }))
}
