## =============================================================================
## longitudinal_heatmap.R
##
## A small toolkit of functions for longitudinal / time-course omics analysis
## (RNA-seq, ATAC-seq, CUT&RUN, drug-treatment time courses, ...).
##
## Given a matrix of features (genes or peaks) measured across an ordered set of
## samples (time points / conditions), these functions cluster features into
## temporal modules, draw ordered heatmaps, summarize module trajectories,
## export module tables, label selected features, and align paired omics layers
## such as RNA with ATAC or CUT&RUN.
##
## The functions are deliberately generic: nothing is hard-coded to a particular
## project.  Feature selection (which genes/peaks are "interesting") is expected
## to be done upstream; you pass in the matrix you want to plot.
##
## Package dependencies are declared in DESCRIPTION. When using this file by
## source(), install the packages listed in README.md first.
## =============================================================================

## Palette registry ------------------------------------------------------------
##
## Two kinds of palette live here and they must be sampled differently:
##
##   - qualitative palettes (Paired, Set2, Set3, Okabe-Ito, ...) are sets of
##     hand-picked, maximally distinguishable colours.  For n <= length(base) we
##     take the FIRST n as-is.  Interpolating them would blend neighbouring hues
##     into muddy in-between colours that are not in the palette at all (and, for
##     Okabe-Ito, would destroy its colour-blind safety).
##   - ramp palettes (Spectral, Viridis) are continuous colour scales, so the
##     correct way to get n colours is to interpolate across the whole ramp.
.moduleviz_palettes <- list(
  paired    = list(ramp = FALSE, cols = NULL),   # filled in below from RColorBrewer
  set2      = list(ramp = FALSE, cols = NULL),
  set3      = list(ramp = FALSE, cols = NULL),
  spectral  = list(ramp = TRUE,  cols = NULL),
  viridis   = list(ramp = TRUE,  cols = NULL),
  okabe_ito = list(ramp = FALSE,
                   cols = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                            "#56B4E9", "#F0E442", "#000000")),
  nature    = list(ramp = FALSE,
                   cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                            "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")),
  science   = list(ramp = FALSE,
                   cols = c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
                            "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"))
)

.moduleviz_palette_names <- function() names(.moduleviz_palettes)

.moduleviz_palette_base <- function(palette) {
  switch(palette,
    paired   = RColorBrewer::brewer.pal(12, "Paired"),
    set2     = RColorBrewer::brewer.pal(8,  "Set2"),
    set3     = RColorBrewer::brewer.pal(12, "Set3"),
    spectral = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    viridis  = grDevices::hcl.colors(11, "viridis"),
    .moduleviz_palettes[[palette]]$cols)
}

#' Build `n` discrete colours from a named ModuleViz palette.
#'
#' @param n number of colours needed.
#' @param palette one of "paired", "set2", "set3", "spectral", "viridis",
#'   "okabe_ito", "nature", "science".
#' @return character vector of `n` hex colours.
#' @noRd
.moduleviz_discrete_palette <- function(n, palette = "paired") {
  palette <- .moduleviz_match_palette(palette)
  base <- .moduleviz_palette_base(palette)
  if (!.moduleviz_palettes[[palette]]$ramp && n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::colorRampPalette(base)(n)
}

## Accept a few friendly aliases (and the historical "spectual" typo).
.moduleviz_match_palette <- function(palette) {
  palette <- tolower(as.character(palette)[1])
  palette <- switch(palette,
    spectual = "spectral", okabeito = "okabe_ito", "okabe-ito" = "okabe_ito",
    rdylbu = "spectral", palette)
  match.arg(palette, .moduleviz_palette_names())
}

.moduleviz_cluster_cols <- function(levels, palette = "paired") {
  setNames(.moduleviz_discrete_palette(length(levels), palette), levels)
}

## Continuous z-score colour maps ---------------------------------------------

#' Diverging / sequential colour map for z-scores.
#'
#' @param palette "rdbu" (default; ColorBrewer RdBu reversed so high = red) or
#'   "viridis" (perceptually uniform, colour-blind safe, good for print).
#' @param limits numeric length-2 z-score range to saturate at.
#' @return a `circlize::colorRamp2` function.
#' @noRd
.moduleviz_zscore_col_fun <- function(palette = c("rdbu", "viridis"),
                                      limits = c(-2, 2)) {
  palette <- tolower(as.character(palette)[1])
  palette <- match.arg(palette, c("rdbu", "viridis"))
  cols <- if (palette == "viridis") {
    grDevices::hcl.colors(7, "viridis")
  } else {
    ## brewer RdBu runs red -> blue; reverse it so LOW = blue, HIGH = red.
    rev(RColorBrewer::brewer.pal(7, "RdBu"))
  }
  circlize::colorRamp2(seq(limits[1], limits[2], length.out = length(cols)), cols)
}

## Back-compat alias for the previous internal name.
.longi_default_col_fun <- function(palette = "rdbu", limits = c(-2, 2)) {
  .moduleviz_zscore_col_fun(palette, limits)
}

## Legend text sizing, shared by every heatmap in the package so that one
## `legend_fontsize` argument moves the title and the labels together.
.moduleviz_legend_gp <- function(legend_fontsize) {
  list(title_gp  = gpar(fontsize = legend_fontsize, fontface = "bold"),
       labels_gp = gpar(fontsize = legend_fontsize))
}

## One ggplot theme for every line / dot plot in the package: thin black panel
## border, no grid, no strip boxes, black text.
.moduleviz_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", size = base_size),
          legend.key.size = unit(4, "mm"),
          legend.background = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          panel.spacing = unit(2.5, "mm"),
          plot.margin = margin(4, 4, 4, 4))
}

## Column labels for a heatmap body.  "time" prints only the time point under
## each column (the group is already shown by the column split / top bar), so
## long "Group|Time" sample names never crowd the figure.
.moduleviz_column_labels <- function(z, meta, time_col, column_labels) {
  labs <- colnames(z)
  is_time <- FALSE
  if (identical(column_labels, "time") && !is.null(meta) &&
      !is.null(time_col) && time_col %in% names(meta) &&
      nrow(meta) == ncol(z)) {
    labs <- as.character(meta[[time_col]])
    is_time <- TRUE
    ## replicates share a time point: print it once per run of columns
    labs[c(FALSE, labs[-1] == labs[-length(labs)])] <- ""
  }
  list(labels = labs, is_time = is_time,
       rot = if (max(nchar(labs)) <= 3) 0 else 90)
}

.moduleviz_write_plot <- function(file, width, height, draw_it) {
  ext <- tolower(sub("^.*\\.([^.]+)$", "\\1", file))
  switch(ext,
    png  = grDevices::png(file, width = width, height = height, units = "in",
                          res = 300),
    jpg  = ,
    jpeg = grDevices::jpeg(file, width = width, height = height, units = "in",
                           res = 300),
    tif  = ,
    tiff = grDevices::tiff(file, width = width, height = height, units = "in",
                           res = 300),
    svg  = grDevices::svg(file, width = width, height = height),
    grDevices::pdf(file, width = width, height = height))
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_it()
  message("wrote ", file)
}

## =============================================================================
## Low-level helpers
## =============================================================================

#' Z-score every row of a matrix (mean 0, sd 1 across that row's samples).
#'
#' Rows with zero variance (or all-NA) become 0 rather than NaN, so they can
#' still be plotted / clustered without blowing up.
#'
#' @param mat numeric matrix, rows = features, columns = samples.
#' @return matrix of the same shape, z-scored by row.
zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  dimnames(z) <- dimnames(mat)
  z
}

#' Average the columns of a matrix within groups (e.g. collapse replicates to a
#' single value per time point / condition).
#'
#' @param mat numeric matrix, columns = samples.
#' @param grp factor/character of length ncol(mat) giving the group of each column.
#' @param group_levels optional ordering of the groups in the output.
#' @return matrix with one column per group, columns ordered by group_levels.
collapse_columns <- function(mat, grp, group_levels = NULL) {
  grp <- as.character(grp)
  if (is.null(group_levels)) group_levels <- unique(grp)
  group_levels <- as.character(group_levels)
  out <- vapply(group_levels, function(g) {
    cols <- which(grp == g)
    if (length(cols) == 0) return(rep(NA_real_, nrow(mat)))
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(mat)))
  matrix(out, nrow = nrow(mat), dimnames = list(rownames(mat), group_levels))
}

## Clustering algorithms offered by longitudinal_cluster() and choose_k().
.moduleviz_methods <- c("kmeans", "hierarchical", "pam", "clara", "fuzzy")

#' Run one clustering, returning an integer cluster label per row.
#'
#' Dispatches to the requested algorithm.  Kept separate from ordering / naming
#' so it can be reused by the stability resampling.
#'
#' @param z numeric matrix (already z-scored), rows = features.
#' @param k number of clusters.
#' @param method one of "kmeans", "hierarchical", "pam".
#' @param seed RNG seed (k-means / PAM are randomised).
#' @param nstart number of random starts for k-means.
#' @param hclust_method linkage for hierarchical clustering (e.g. "ward.D2").
#' @param dist_method distance for hierarchical / PAM ("euclidean", "correlation", ...).
#' @return integer vector of length nrow(z).
.cluster_once <- function(z, k, method = "kmeans", seed = 1,
                          nstart = 50, hclust_method = "ward.D2",
                          dist_method = "euclidean") {
  method <- match.arg(method, .moduleviz_methods)
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1)
    stop("`k` must be a single positive integer.", call. = FALSE)
  n_uniq <- nrow(unique(z))
  k <- min(as.integer(k), n_uniq)
  if (k < 1) stop("Need at least one distinct row to cluster.")
  ## a single module needs no algorithm (and kmeans/pam reject k = n)
  if (k == 1L) return(rep(1L, nrow(z)))
  if (k == nrow(z)) return(seq_len(nrow(z)))

  ## "correlation" distance = 1 - Pearson correlation between feature profiles.
  ## A constant profile has no defined correlation; treat it as uncorrelated
  ## with everything else rather than letting NA reach hclust()/pam().
  make_dist <- function(mat) {
    if (dist_method == "correlation") {
      r <- suppressWarnings(cor(t(mat)))
      r[is.na(r)] <- 0
      diag(r) <- 1
      as.dist(1 - r)
    } else {
      dist(mat, method = dist_method)
    }
  }

  if (method == "kmeans") {
    ## Hartigan-Wong (R's default) rarely leaves empty clusters, which makes it
    ## more robust than Lloyd on tightly structured time-course data.
    km <- .moduleviz_with_seed(seed,
      kmeans(z, centers = k, nstart = nstart, iter.max = 1000))
    return(as.integer(km$cluster))
  }
  if (method == "hierarchical") {
    hc <- hclust(make_dist(z), method = hclust_method)
    return(as.integer(cutree(hc, k = k)))
  }
  if (method == "pam") {
    ## partitioning around medoids -- robust k-means alternative
    pm <- .moduleviz_with_seed(seed,
      cluster::pam(make_dist(z), k = k, diss = TRUE))
    return(as.integer(pm$clustering))
  }
  if (method == "clara") {
    ## PAM on repeated sub-samples: scales to tens of thousands of features,
    ## where pam()/hclust() would need an n x n distance matrix.
    metric <- if (identical(dist_method, "manhattan")) "manhattan" else "euclidean"
    cl <- .moduleviz_with_seed(seed,
      cluster::clara(z, k = k, metric = metric, samples = 50,
                     sampsize = min(nrow(z), 40 + 2 * k), pamLike = TRUE,
                     rngR = TRUE))
    return(as.integer(cl$clustering))
  }
  ## fuzzy c-means (cluster::fanny): every feature gets a membership in every
  ## module; the hard label is the module with the highest membership, which
  ## is returned alongside as the "membership" attribute.
  if (k > nrow(z) / 2 - 1)
    stop("method = \"fuzzy\" needs k <= n/2 - 1 (n = ", nrow(z),
         " features); lower `k`.", call. = FALSE)
  fz <- suppressWarnings(
    cluster::fanny(make_dist(z), k = k, diss = TRUE, memb.exp = 1.5,
                   maxit = 2000))
  structure(as.integer(fz$clustering),
            membership = apply(fz$membership, 1, max))
}

#' Mean profile of every cluster: a clusters x columns matrix.
#'
#' Always returns a matrix, including when there is a single cluster or a
#' single column (where sapply()/t() would silently change orientation).
#' @noRd
.moduleviz_centroids <- function(mat, cluster, labs) {
  out <- vapply(labs, function(cl)
    colMeans(mat[cluster == cl, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(mat)))
  out <- matrix(out, nrow = ncol(mat))
  out[is.nan(out)] <- NA_real_
  out <- t(out)
  dimnames(out) <- list(as.character(labs), colnames(mat))
  out
}

#' Order clusters so the heatmap reads from "high early" to "high late".
#'
#' For each cluster we take its centroid across the ordered time points, find
#' the time point where it peaks (arg-max), and rank clusters by that peak time.
#' Ties are broken by the overall slope (end minus start) so that, among modules
#' peaking at the same time, rising ones come before falling ones.
#'
#' @param z numeric z-scored matrix.
#' @param cluster integer cluster labels (one per row of z).
#' @param time_group factor/character length ncol(z): the time point of each column.
#' @param time_levels ordered unique time points (first = earliest).
#' @return integer vector: the new rank (1 = peaks earliest) for each input label.
.order_clusters_temporal <- function(z, cluster, time_group, time_levels) {
  tp_mean  <- collapse_columns(z, time_group, time_levels)     # feature x timepoint
  labs     <- sort(unique(cluster))
  centers  <- .moduleviz_centroids(tp_mean, cluster, labs)

  ## which.max() on an all-NA row returns integer(0); rank those modules last
  peak_at <- vapply(seq_len(nrow(centers)), function(i) {
    w <- which.max(centers[i, ])
    if (length(w)) w else ncol(centers) + 1L
  }, integer(1))
  slope   <- centers[, ncol(centers)] - centers[, 1]
  slope[is.na(slope)] <- 0
  ord     <- order(peak_at, slope)

  new_rank <- integer(length(labs))
  new_rank[ord] <- seq_along(labs)
  names(new_rank) <- labs
  new_rank[as.character(cluster)]
}

## =============================================================================
## Clustering stability
## =============================================================================

#' Assess how stable each cluster is under row resampling.
#'
#' Repeatedly sub-samples the features, re-clusters, and measures how well each
#' original cluster is recovered (mean maximum Jaccard similarity, in the spirit
#' of fpc::clusterboot).  Values near 1 mean the module is highly reproducible;
#' values below ~0.6 mean it is unstable and should be interpreted with caution.
#'
#' @param z z-scored matrix.
#' @param cluster the reference integer cluster labels on the full data.
#' @param k,method,nstart,hclust_method,dist_method clustering settings (must
#'   match those used for the reference clustering).
#' @param n_resample number of resampling iterations.
#' @param subsample_frac fraction of rows kept in each iteration.
#' @param seed base RNG seed.
#' @return data.table with columns Cluster, MeanJaccard, N (cluster size).
cluster_stability <- function(z, cluster, k, method = "kmeans",
                              nstart = 50, hclust_method = "ward.D2",
                              dist_method = "euclidean",
                              n_resample = 50, subsample_frac = 0.8,
                              seed = 1) {
  if (length(cluster) != nrow(z))
    stop("`cluster` must have one label per row of `z`.", call. = FALSE)
  if (subsample_frac <= 0 || subsample_frac > 1)
    stop("`subsample_frac` must be in (0, 1].", call. = FALSE)
  if (is.factor(cluster)) cluster <- as.character(cluster)
  n     <- nrow(z)
  labs  <- sort(unique(cluster))
  jacc  <- matrix(NA_real_, nrow = length(labs), ncol = n_resample,
                  dimnames = list(as.character(labs), NULL))

  for (b in seq_len(n_resample)) {
    iter_seed <- if (is.null(seed)) NULL else seed + b
    idx <- .moduleviz_with_seed(iter_seed,
      sample.int(n, size = min(n, max(2L, floor(subsample_frac * n)))))
    zb  <- z[idx, , drop = FALSE]
    ## reclustering the subsample; use a per-iteration seed so starts differ.
    ## Empty-cluster warnings are expected on small resamples and handled by
    ## kmeans internally, so silence them here to keep the console readable.
    cb  <- tryCatch(
      suppressWarnings(.cluster_once(zb, k = k, method = method, seed = iter_seed,
                    nstart = nstart, hclust_method = hclust_method,
                    dist_method = dist_method)),
      error = function(e) rep(1L, nrow(zb)))
    ref <- cluster[idx]

    for (li in seq_along(labs)) {
      in_ref <- ref == labs[li]
      if (!any(in_ref)) next
      ## best Jaccard overlap with any resample cluster
      jacc[li, b] <- max(vapply(sort(unique(cb)), function(d) {
        in_new <- cb == d
        inter  <- sum(in_ref & in_new)
        union  <- sum(in_ref | in_new)
        if (union == 0) 0 else inter / union
      }, numeric(1)))
    }
  }

  data.table(
    Cluster     = if (is.numeric(labs)) paste0("Cluster", labs) else as.character(labs),
    MeanJaccard = round(rowMeans(jacc, na.rm = TRUE), 3),
    N           = as.integer(table(factor(cluster, levels = labs)))
  )
}

#' Scan a range of k and report elbow / silhouette diagnostics.
#'
#' Helps choose the number of clusters before committing.  Returns a table and,
#' optionally, saves a diagnostic plot.  The "elbow" is the k furthest from the
#' straight line joining the first and last (k, within-ss) points.
#'
#' @param z z-scored matrix.
#' @param k_range integer vector of k values to try.
#' @param method clustering method (silhouette is computed for all methods).
#' @param seed base seed.
#' @param file optional path; if given, a PDF diagnostic plot is written.
#' @param chosen_k optional k actually used downstream.  When supplied it is
#'   drawn as a solid line beside the dashed elbow, so the diagnostic figure
#'   stays consistent with the module count in the rest of the analysis - the
#'   elbow is only a heuristic and is often not the k you want to report.
#' @param annotate_k if FALSE, draw the bare elbow curve with no marker lines
#'   and no subtitle.
#' @param base_size base font size for the diagnostic plot.
#' @return data.table of diagnostics; attribute "suggested_k" holds the elbow k.
choose_k <- function(z, k_range = 2:20, method = "kmeans",
                     nstart = 25, hclust_method = "ward.D2",
                     dist_method = "euclidean", seed = 100, file = NULL,
                     chosen_k = NULL, annotate_k = TRUE, base_size = 10) {
  z <- as.matrix(z)
  k_range <- sort(unique(as.integer(k_range)))
  k_range <- k_range[k_range >= 1 & k_range < nrow(unique(z))]
  if (!length(k_range))
    stop("No usable k in `k_range`: every value must be smaller than the ",
         "number of distinct rows (", nrow(unique(z)), ").", call. = FALSE)
  d_for_sil <- if (dist_method == "correlation") {
    r <- suppressWarnings(cor(t(z)))
    r[is.na(r)] <- 0
    diag(r) <- 1
    as.dist(1 - r)
  } else {
    dist(z, method = dist_method)
  }

  diag_dt <- rbindlist(lapply(k_range, function(k) {
    cl <- .cluster_once(z, k = k, method = method,
                        seed = if (is.null(seed)) NULL else seed + k,
                        nstart = nstart, hclust_method = hclust_method,
                        dist_method = dist_method)
    ## within-cluster sum of squares (to the cluster mean)
    wss <- sum(vapply(sort(unique(cl)), function(c) {
      sub <- z[cl == c, , drop = FALSE]
      if (nrow(sub) < 2) return(0)
      sum(scale(sub, center = TRUE, scale = FALSE)^2)
    }, numeric(1)))
    sil <- NA_real_
    if (length(unique(cl)) > 1 && length(unique(cl)) < nrow(z)) {
      sil <- mean(cluster::silhouette(cl, d_for_sil)[, "sil_width"])
    }
    data.table(k = k, tot_withinss = wss, mean_silhouette = sil)
  }))

  ## distance-to-line elbow
  x <- diag_dt$k; y <- diag_dt$tot_withinss
  x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  ## with one or two candidate k there is no curve to find an elbow on
  elbow <- if (length(x) < 3 || !is.finite(denom) || denom == 0) {
    x[1]
  } else {
    dl <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / denom
    x[which.max(dl)]
  }
  diag_dt[, elbow_k := elbow]
  setattr(diag_dt, "suggested_k", elbow)

  if (!is.null(file)) {
    subtitle <- NULL
    p <- ggplot(diag_dt, aes(k, tot_withinss))
    if (annotate_k) {
      subtitle <- sprintf("dashed: elbow (k = %s)", elbow)
      p <- p + geom_vline(xintercept = elbow, linetype = 2, color = "#b2182b")
      if (!is.null(chosen_k)) {
        p <- p + geom_vline(xintercept = chosen_k, linetype = 1,
                            linewidth = 0.5, color = "#2166ac")
        subtitle <- sprintf("%s; solid: chosen (k = %s)", subtitle, chosen_k)
      }
    }
    p <- p +
      geom_line(color = "black", linewidth = 0.5) +
      geom_point(size = 1.5, color = "black") +
      scale_x_continuous(breaks = k_range) +
      .moduleviz_theme(base_size) +
      theme(plot.subtitle = element_text(size = base_size * 0.85)) +
      labs(x = "Number of modules (k)", y = "Within-cluster sum of squares",
           subtitle = subtitle)
    ggsave(file, p, width = 4, height = 3)
  }
  diag_dt[]
}

## =============================================================================
## Longitudinal clustering
## =============================================================================

#' Cluster longitudinal features into temporal modules.
#'
#' This is the workhorse.  It z-scores the matrix (unless told the input is
#' already z-scored), clusters the features, orders the clusters from
#' "high early" to "high late", orders the rows within each cluster by signal
#' strength, and (optionally) assesses cluster stability.
#'
#' @param mat matrix / data.frame / data.table of features x samples.  If a
#'   data.table/data.frame with a character ID column is passed, set `id_col`.
#' @param meta data.frame with one row per sample.  Must contain `sample_col`
#'   and `time_col`; may contain `group_col` (e.g. genotype / treatment arm).
#' @param sample_col column in `meta` giving sample IDs (matching mat columns).
#' @param time_col column in `meta` giving a *numeric* time (used for ordering).
#' @param group_col optional column in `meta` for a second factor (e.g. genotype)
#'   that will split the heatmap columns and the line-plot facets.
#' @param id_col if `mat` carries feature IDs in a column, its name.
#' @param k number of clusters/modules.
#' @param method "kmeans" (default), "hierarchical", or "pam".
#' @param is_zscore set TRUE if `mat` is already z-scored (skip scaling).
#' @param aggregate_replicates if TRUE, average replicate columns within each
#'   group_col x time_col before clustering (cleaner, fewer columns).
#' @param seed,nstart,hclust_method,dist_method clustering controls.
#' @param stability if TRUE, run cluster_stability().
#' @param n_resample,subsample_frac stability controls.
#' @return object of class "longi": a list with z (ordered matrix), cluster
#'   (ordered factor), membership, centroids, stability, meta (ordered), params.
longitudinal_cluster <- function(mat, meta = NULL,
                                 sample_col = "SampleID",
                                 time_col   = "TimeHr",
                                 group_col  = NULL,
                                 id_col     = NULL,
                                 k = 8,
                                 method = "kmeans",
                                 is_zscore = FALSE,
                                 aggregate_replicates = FALSE,
                                 seed = 1, nstart = 50,
                                 hclust_method = "ward.D2",
                                 dist_method = "euclidean",
                                 stability = TRUE,
                                 n_resample = 50, subsample_frac = 0.8,
                                 assay = NULL) {

  method <- match.arg(method, .moduleviz_methods)

  ## ---- accept a SummarizedExperiment / SingleCellExperiment ---------------
  input <- .moduleviz_unpack(mat, meta, sample_col = sample_col, assay = assay)
  mat <- input$mat
  meta <- input$meta
  sample_col <- input$sample_col
  if (is.null(meta))
    stop("`meta` (one row per sample) is required unless `mat` is a ",
         "SummarizedExperiment.", call. = FALSE)

  meta <- as.data.table(copy(as.data.frame(meta, stringsAsFactors = FALSE)))
  need <- c(sample_col, time_col, group_col)
  absent <- setdiff(need, names(meta))
  if (length(absent))
    stop("Column(s) not found in `meta`: ", paste(absent, collapse = ", "),
         ". Available: ", paste(names(meta), collapse = ", "), call. = FALSE)
  if (anyDuplicated(as.character(meta[[sample_col]])))
    stop("`meta[[\"", sample_col, "\"]]` contains duplicated sample IDs.",
         call. = FALSE)

  ## ---- build the numeric matrix with feature rownames --------------------
  sample_ids <- as.character(meta[[sample_col]])
  if (!is.null(id_col)) {
    dt  <- as.data.table(mat)
    if (!id_col %in% names(dt))
      stop("`id_col` '", id_col, "' is not a column of `mat`.", call. = FALSE)
    ids <- as.character(dt[[id_col]])
    keep_cols <- intersect(sample_ids, names(dt))
    m   <- as.matrix(dt[, keep_cols, with = FALSE])
    rownames(m) <- ids
  } else {
    m <- as.matrix(mat)
    if (is.null(rownames(m)))
      stop("mat has no rownames; supply feature IDs as rownames or via id_col.")
    keep_cols <- intersect(sample_ids, colnames(m))
    m <- m[, keep_cols, drop = FALSE]
  }
  if (length(keep_cols) == 0)
    stop("No overlap between meta[[sample_col]] and the matrix columns.")
  if (!is.numeric(m))
    stop("`mat` must be numeric; non-numeric columns were found. Use `id_col` ",
         "if feature IDs are stored in a column.", call. = FALSE)
  if (anyDuplicated(rownames(m))) {
    warning("Duplicated feature IDs were made unique with make.unique().",
            call. = FALSE)
    rownames(m) <- make.unique(rownames(m))
  }

  ## keep meta rows matching present samples, in matrix-column order
  meta <- meta[match(keep_cols, sample_ids)]
  if (anyNA(meta[[time_col]]))
    stop("`meta[[\"", time_col, "\"]]` contains missing values.", call. = FALSE)
  complete <- complete.cases(m) & rowSums(!is.finite(m)) == 0
  if (!all(complete))
    message("Dropping ", sum(!complete),
            " feature(s) with missing or non-finite values.")
  m    <- m[complete, , drop = FALSE]
  if (nrow(m) < 2) stop("Fewer than 2 complete-case features to cluster.")

  ## ---- order samples: group (if any) then time then remaining ------------
  ## Time may be numeric, a factor (level order), or character (order of first
  ## appearance); `time_levels` is that order as character keys.
  time_levels <- .moduleviz_time_levels(meta[[time_col]])
  time_rank   <- match(as.character(meta[[time_col]]), time_levels)
  ord_keys <- if (is.null(group_col)) list(time_rank)
              else list(meta[[group_col]], time_rank)
  col_order <- do.call(order, c(ord_keys, list(seq_len(nrow(meta)))))
  meta <- meta[col_order]
  m    <- m[, as.character(meta[[sample_col]]), drop = FALSE]

  ## optional replicate collapsing (one column per group x time) ------------
  if (aggregate_replicates) {
    key <- if (is.null(group_col)) as.character(meta[[time_col]])
           else paste(meta[[group_col]], meta[[time_col]], sep = "|")
    key_levels <- unique(key)
    m <- collapse_columns(m, key, key_levels)
    ## One metadata row per aggregated column.  Group and time keep their
    ## original type, and any other column that is constant within a
    ## group x time cell (e.g. a condition encoded in the group) is carried
    ## over; per-replicate columns cannot be and are dropped.
    carry <- setdiff(.moduleviz_constant_cols(meta, key), sample_col)
    agg_meta <- meta[match(key_levels, key), carry, with = FALSE]
    set(agg_meta, j = sample_col, value = key_levels)
    setcolorder(agg_meta, sample_col)
    meta <- agg_meta
  }

  ## ---- z-score ------------------------------------------------------------
  z <- if (is_zscore) m else zscore_rows(m)

  ## ---- cluster ------------------------------------------------------------
  raw <- .cluster_once(z, k = k, method = method, seed = seed,
                       nstart = nstart, hclust_method = hclust_method,
                       dist_method = dist_method)
  k_use <- length(unique(raw))

  ## time grouping for ordering: the time point of each column
  time_group  <- as.character(meta[[time_col]])
  rank_vec    <- .order_clusters_temporal(z, raw, time_group, time_levels)

  cluster <- factor(rank_vec,
                    levels = seq_len(k_use),
                    labels = paste0("Cluster", seq_len(k_use)))

  ## ---- order rows: by cluster, then by strongest signal ------------------
  row_score <- apply(abs(z), 1, max)
  row_order <- order(cluster, -row_score)
  z         <- z[row_order, , drop = FALSE]
  cluster   <- cluster[row_order]
  row_score <- row_score[row_order]

  membership <- data.table(ID = rownames(z),
                           Cluster = as.character(cluster),
                           MaxAbsZ = round(row_score, 3))
  ## fuzzy clustering also reports how firmly each feature sits in its module
  if (!is.null(attr(raw, "membership")))
    membership[, Membership := round(attr(raw, "membership")[row_order], 3)]

  ## ---- cluster centroids per (group x) time (for the line plot) ----------
  if (is.null(group_col)) {
    tp_mean  <- collapse_columns(z, as.character(meta[[time_col]]), time_levels)
  } else {
    gkey <- paste(meta[[group_col]], meta[[time_col]], sep = "|")
    tp_mean <- collapse_columns(z, gkey, unique(gkey))
  }
  centroid <- .moduleviz_centroids(tp_mean, cluster, levels(cluster))
  centroid_dt <- data.table(Cluster = rownames(centroid), as.data.table(centroid))

  ## ---- stability ----------------------------------------------------------
  stab <- NULL
  if (stability) {
    ## stability uses the (unordered) raw labels on the *ordered* z; recompute
    ## raw labels aligned to ordered z from the membership so mapping is exact.
    raw_ordered <- as.integer(sub("Cluster", "", as.character(cluster)))
    stab <- cluster_stability(z, raw_ordered, k = k_use, method = method,
                              nstart = nstart, hclust_method = hclust_method,
                              dist_method = dist_method,
                              n_resample = n_resample,
                              subsample_frac = subsample_frac, seed = seed)
  }

  structure(list(
    z          = z,
    cluster    = cluster,
    membership = membership,
    centroids  = centroid_dt,
    stability  = stab,
    meta       = meta,
    params     = list(k = k, k_use = k_use, method = method,
                      sample_col = sample_col, time_col = time_col,
                      group_col = group_col, seed = seed,
                      dist_method = dist_method, hclust_method = hclust_method,
                      time_levels = time_levels)
  ), class = "longi")
}

print.longi <- function(x, ...) {
  p <- x$params
  cat(sprintf("<longi>  %d features x %d samples,  %d modules (%s)\n",
              nrow(x$z), ncol(x$z), p$k_use, p$method))
  if (!is.null(x$stability))
    cat(sprintf("  stability (mean Jaccard): %.2f - %.2f\n",
                min(x$stability$MeanJaccard), max(x$stability$MeanJaccard)))
  invisible(x)
}

## =============================================================================
## 5. Choosing which features to label (helper used by the heatmap)
## =============================================================================

#' Resolve the set of features to label on a heatmap into row indices.
#'
#' Supports three ways to choose labels, combined in this order:
#'   - `features`: an explicit character vector.  Each entry is matched first as
#'     an exact ID, then (if not found) as a gene-symbol / substring match.
#'   - `top_n`: automatically take the top-N features per cluster, ranked by
#'     `score` (default = max |z|, i.e. strongest signal).
#'   - `rank_table`: an optional data.table of (ID, value) to rank by instead
#'     (e.g. -log10 padj or |log2FC| from a DE test).
#'
#' @param obj a "longi" object.
#' @param features explicit IDs/symbols to label (character vector) or NULL.
#' @param top_n integer; per-cluster automatic labels (0 = none).
#' @param rank_table optional data.table with columns id_field & value_field.
#' @param id_field,value_field column names in rank_table.
#' @param bigger_is_better rank direction for rank_table / score.
#' @return data.table(idx, label) of row positions and label text.
resolve_labels <- function(obj, features = NULL, top_n = 0,
                           rank_table = NULL, id_field = "ID",
                           value_field = "value", bigger_is_better = TRUE) {
  ids <- rownames(obj$z)
  out <- list()

  ## (a) explicit, user-chosen features
  if (!is.null(features) && length(features)) {
    idx <- integer(0); lab <- character(0)
    ## Match in order of decreasing strictness: the exact ID, then a whole
    ## "_"/"|"-delimited token (so "Il2" finds "peak7_Il2" but not
    ## "peak9_Il2ra"), and only then a plain substring.
    tokens <- strsplit(ids, "[_|]")
    unmatched <- character(0)
    for (f in as.character(features)) {
      hit <- which(ids == f)
      if (!length(hit))
        hit <- which(vapply(tokens, function(tk) f %in% tk, logical(1)))
      if (!length(hit)) hit <- grep(f, ids, fixed = TRUE)
      if (length(hit)) { idx <- c(idx, hit); lab <- c(lab, ids[hit]) }
      else unmatched <- c(unmatched, f)
    }
    if (length(unmatched))
      message("Features to label not found: ",
              paste(unmatched, collapse = ", "))
    if (length(idx)) out[[length(out) + 1]] <- data.table(idx = idx, label = lab)
  }

  ## (b) automatic top-N per cluster
  if (top_n > 0) {
    dt <- data.table(idx = seq_along(ids), ID = ids,
                     Cluster = as.character(obj$cluster))
    if (!is.null(rank_table)) {
      rt <- as.data.table(rank_table)[, .(ID = get(id_field), value = get(value_field))]
      dt <- merge(dt, rt, by = "ID", all.x = TRUE)
    } else {
      dt[, value := apply(abs(obj$z), 1, max)]
    }
    setorderv(dt, "value", if (bigger_is_better) -1L else 1L)
    dt <- dt[!is.na(value)]
    picked <- dt[, head(.SD, top_n), by = Cluster]
    out[[length(out) + 1]] <- picked[, .(idx, label = ID)]
  }

  if (!length(out)) return(data.table(idx = integer(0), label = character(0)))
  res <- unique(rbindlist(out))
  setorder(res, idx)
  res
}

## =============================================================================
## Ordered heatmap
## =============================================================================

#' Draw the ordered longitudinal heatmap.
#'
#' Rows are split by module (ordered early-peaking -> late-peaking); columns are
#' split by group_col (if set) and ordered by time.  Metadata columns become
#' coloured annotation bars along the top.  Chosen features are labelled on the
#' right via link lines.
#'
#' @param obj a "longi" object from longitudinal_cluster().
#' @param annotation_cols character vector of meta columns to show as top bars.
#'   Defaults to group_col, plus time_col when the time points are not already
#'   printed under the columns.
#' @param annotation_colors optional named list mapping annotation -> colour vec.
#' @param label_features explicit feature IDs/symbols to label (see resolve_labels).
#' @param top_n_label per-cluster automatic labels (0 = none).
#' @param rank_table optional ranking table for automatic labels.
#' @param id_field,value_field,bigger_is_better passed to resolve_labels.
#' @param col_fun colour mapping for z-scores (circlize::colorRamp2).  Overrides
#'   `zscore_palette` when supplied.
#' @param zscore_palette continuous palette for the z-score body: "rdbu"
#'   (default, blue-white-red) or "viridis".
#' @param zscore_limits numeric length-2 z-score range to saturate the colour
#'   scale at.
#' @param cluster_palette discrete palette for the module annotation bar.
#'   Defaults to "paired"; "spectral" and the other registry palettes also work.
#' @param annotation_palette discrete palette for categorical metadata bars
#'   (e.g. condition / genotype).  Defaults to "set2"; "set3" is a good choice
#'   when there are more than 8 levels.
#' @param label_fontsize,row_title_fontsize,column_name_fontsize text sizes.
#' @param annotation_name_fontsize,legend_fontsize,title_fontsize text sizes for
#'   the annotation track names, all legends, and the heatmap title.
#' @param show_column_names show sample names along the bottom.  Set FALSE when
#'   sample IDs are long or uninformative and the annotation bars already say
#'   what each column is.
#' @param column_labels what to print under each column: "time" (default) the
#'   time point only, since the group is already shown by the column split and
#'   the top bar; "sample" the full sample names.
#' @param show_stability if TRUE and available, add a module stability bar
#'   (off by default; the values are in the stability table).
#' @param file optional PDF or PNG path; if set the heatmap is written there.
#' @param width,height PDF size in inches.
#' @param title optional heatmap title.
#' @return the (invisibly returned) Heatmap object; drawn to the device / file.
longitudinal_heatmap <- function(obj,
                                 annotation_cols = NULL,
                                 annotation_colors = NULL,
                                 label_features = NULL,
                                 top_n_label = 0,
                                 rank_table = NULL,
                                 id_field = "ID", value_field = "value",
                                 bigger_is_better = TRUE,
                                 col_fun = NULL,
                                 zscore_palette = c("rdbu", "viridis"),
                                 zscore_limits = c(-2, 2),
                                 cluster_palette = "paired",
                                 annotation_palette = "set2",
                                 label_fontsize = 9,
                                 row_title_fontsize = 10,
                                 column_name_fontsize = 9,
                                 annotation_name_fontsize = 9,
                                 legend_fontsize = 9,
                                 title_fontsize = 11,
                                 show_column_names = TRUE,
                                 show_stability = FALSE,
                                 file = NULL, width = 5, height = 6.5,
                                 title = NULL,
                                 enrichment = NULL,
                                 enrichment_top_n = 3,
                                 enrichment_padj = 0.05,
                                 enrichment_fontsize = 9,
                                 column_labels = c("time", "sample")) {
  stopifnot(inherits(obj, "longi"))
  column_labels <- match.arg(column_labels)
  if (is.null(col_fun))
    col_fun <- .moduleviz_zscore_col_fun(zscore_palette, zscore_limits)
  p    <- obj$params
  meta <- obj$meta
  z    <- obj$z
  cluster <- obj$cluster

  ## ---- column split by group (if provided) --------------------------------
  if (!is.null(p$group_col)) {
    col_split <- factor(meta[[p$group_col]], levels = unique(meta[[p$group_col]]))
  } else {
    col_split <- NULL
  }

  ## ---- top annotation bars ------------------------------------------------
  clab <- .moduleviz_column_labels(z, meta, p$time_col, column_labels)
  if (is.null(annotation_cols)) {
    ## time is on the axis already, so it does not need a bar of its own
    annotation_cols <- c(p$group_col,
                         if (!(clab$is_time && show_column_names)) p$time_col)
    annotation_cols <- annotation_cols[!vapply(annotation_cols, is.null, logical(1))]
  }
  ## Requested-but-absent annotations are dropped, but say so: after
  ## aggregate_replicates = TRUE the metadata is rebuilt with only the sample,
  ## group, and time columns, so any other column silently disappears.
  missing_anno <- setdiff(annotation_cols, names(meta))
  if (length(missing_anno))
    warning("annotation_cols not found in the object's metadata and skipped: ",
            paste(missing_anno, collapse = ", "),
            ". Available: ", paste(names(meta), collapse = ", "), call. = FALSE)
  annotation_cols <- intersect(annotation_cols, names(meta))
  top_anno <- NULL
  if (length(annotation_cols)) {
    anno_df <- as.data.frame(meta[, annotation_cols, with = FALSE])
    ## build colours for any annotation the user didn't specify
    col_list <- if (is.null(annotation_colors)) list() else annotation_colors
    for (a in annotation_cols) {
      if (!is.null(col_list[[a]])) next
      vals <- anno_df[[a]]
      if (is.numeric(vals)) {
        ## sequential ramp for continuous metadata (e.g. time)
        col_list[[a]] <- circlize::colorRamp2(
          range(vals, na.rm = TRUE),
          c("#deebf7", "#08519c"))
      } else {
        lv <- unique(as.character(vals))
        col_list[[a]] <- setNames(
          .moduleviz_discrete_palette(length(lv), annotation_palette), lv)
      }
    }
    top_anno <- HeatmapAnnotation(
      df = anno_df, col = col_list,
      simple_anno_size = unit(3, "mm"),
      annotation_name_gp = gpar(fontsize = annotation_name_fontsize),
      annotation_legend_param = .moduleviz_legend_gp(legend_fontsize))
  }

  ## ---- left annotation: module colour + optional stability ---------------
  cluster_cols <- .moduleviz_cluster_cols(levels(cluster), cluster_palette)
  ## the row titles already name the modules, so the bar needs no legend
  left_args <- list(Cluster = cluster,
                    col = list(Cluster = cluster_cols),
                    show_annotation_name = FALSE,
                    show_legend = FALSE,
                    annotation_legend_param = .moduleviz_legend_gp(legend_fontsize),
                    simple_anno_size = unit(2.5, "mm"))
  if (show_stability && !is.null(obj$stability)) {
    left_args$show_legend <- c(FALSE, TRUE)
    stab_map <- setNames(obj$stability$MeanJaccard, obj$stability$Cluster)
    left_args$Stability <- stab_map[as.character(cluster)]
    left_args$col$Stability <- circlize::colorRamp2(c(0.5, 0.75, 1),
                                                    c("#f7f7f7", "#bdbdbd", "#252525"))
  }
  left_anno <- do.call(rowAnnotation, left_args)

  ## ---- row titles with module size ----------------------------------------
  cn <- table(cluster)
  row_title <- sprintf("%s\n(n=%d)", levels(cluster), as.integer(cn))

  ## ---- right annotation: feature labels -----------------------------------
  lab <- resolve_labels(obj, features = label_features, top_n = top_n_label,
                        rank_table = rank_table, id_field = id_field,
                        value_field = value_field,
                        bigger_is_better = bigger_is_better)
  right_parts <- list()
  if (nrow(lab)) {
    right_parts$mark <- anno_mark(which = "row",
      at = lab$idx, labels = lab$label,
      labels_gp = gpar(fontsize = label_fontsize, fontface = "italic"),
      link_gp = gpar(lwd = 0.5), link_width = unit(4, "mm"))
  }
  ## ---- right annotation: top enriched terms per module --------------------
  if (!is.null(enrichment)) {
    right_parts$Enrichment <- .moduleviz_enrichment_anno(
      enrichment, levels(cluster), top_n = enrichment_top_n,
      padj_cutoff = enrichment_padj, fontsize = enrichment_fontsize)
  }
  right_anno <- NULL
  if (length(right_parts))
    right_anno <- do.call(rowAnnotation,
                          c(right_parts, list(show_annotation_name = FALSE)))

  ht <- Heatmap(
    z, name = "z-score", col = col_fun,
    row_split = cluster, cluster_rows = FALSE, cluster_row_slices = FALSE,
    row_title = row_title, row_title_rot = 0,
    row_title_gp = gpar(fontsize = row_title_fontsize),
    row_gap = unit(1, "mm"), show_row_names = FALSE,
    column_split = col_split, cluster_columns = FALSE,
    column_gap = unit(1, "mm"),
    show_column_names = show_column_names,
    column_labels = clab$labels, column_names_rot = clab$rot,
    column_names_centered = clab$rot == 0,
    column_names_gp = gpar(fontsize = column_name_fontsize),
    ## the 6 cm default clips long sample names at larger font sizes
    column_names_max_height = max_text_width(
      clab$labels, gp = gpar(fontsize = column_name_fontsize)),
    column_title = title,
    column_title_gp = gpar(fontsize = title_fontsize),
    top_annotation = top_anno, left_annotation = left_anno,
    right_annotation = right_anno,
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = c(list(direction = "horizontal",
                                  title_position = "topcenter",
                                  legend_width = unit(2.5, "cm"),
                                  grid_height = unit(2.5, "mm")),
                             .moduleviz_legend_gp(legend_fontsize)))

  draw_it <- function() draw(ht, heatmap_legend_side = "bottom",
                             annotation_legend_side = "bottom",
                             merge_legend = TRUE)
  if (!is.null(file)) {
    .moduleviz_write_plot(file, width, height, draw_it)
  } else {
    draw_it()
  }
  invisible(ht)
}

## =============================================================================
## Module trajectory plot
## =============================================================================

#' Line plot of the mean z-score trajectory per module.
#'
#' One facet per module (rows).  If the "longi" object has a group_col, the
#' groups are drawn as separate coloured lines / facet columns so you can
#' compare, e.g., WT vs KO trajectories side by side.
#'
#' @param obj a "longi" object.
#' @param file optional PDF path.
#' @param width,height PDF size; height auto-scales with the module count if NULL.
#' @param ncol number of facet columns when there is no group_col.
#' @param palette colour palette for the group / condition lines.  Defaults to
#'   "set2"; any registry palette works ("set3", "paired", "nature", ...).
#' @param cluster_palette palette used to colour the per-module lines when the
#'   object has no group_col.  Defaults to "paired".
#' @param group_colors optional named character vector of explicit colours for
#'   the group lines, e.g. `c(Control = "grey70", Stimulated = "firebrick")`.
#'   Overrides `palette` for the levels it names; any level it does not name
#'   keeps its palette colour.
#' @param base_size ggplot base font size.
#' @param title optional plot title; none by default.
#' @param title_fontsize title size in points.  Defaults to `base_size`, so the
#'   title stays proportionate instead of inflating with the base size.
#' @param xlab x-axis label; defaults to the name of `time_col`.
#' @param time_labels optional relabelling of the x-axis ticks, e.g.
#'   `c("1" = "Naive", "2" = "Memory", "3" = "Exhausted")` when the ordering
#'   column is an integer rank rather than real time.
#' @param x_angle rotation of the x tick labels in degrees.  Use `0` when the
#'   labels are short enough to sit flat.
#' @param line_width,point_size line and point size.
#' @return a ggplot object (also saved if `file` is given).
pattern_lineplot <- function(obj, file = NULL, width = NULL, height = NULL,
                             ncol = 1,
                             palette = "set2",
                             cluster_palette = "paired",
                             group_colors = NULL,
                             base_size = 10,
                             title = NULL,
                             title_fontsize = NULL,
                             xlab = NULL,
                             time_labels = NULL,
                             x_angle = 0,
                             line_width = 0.7,
                             point_size = 1.5) {
  stopifnot(inherits(obj, "longi"))
  palette <- .moduleviz_match_palette(palette)
  cluster_palette <- .moduleviz_match_palette(cluster_palette)
  if (is.null(title_fontsize)) title_fontsize <- base_size
  p  <- obj$params
  cn <- table(obj$cluster)
  facet_lab <- setNames(sprintf("%s (n=%d)", levels(obj$cluster), as.integer(cn)),
                        levels(obj$cluster))

  long <- melt(obj$centroids, id.vars = "Cluster",
               variable.name = "key", value.name = "z")
  long[, Cluster := factor(Cluster, levels = levels(obj$cluster))]
  long[, ClusterF := factor(facet_lab[as.character(Cluster)],
                            levels = facet_lab[levels(obj$cluster)])]

  ## Recover group and time for every centroid column from the metadata
  ## rather than by parsing the "Group|Time" keys, so non-numeric time points
  ## and group names containing "|" both work.
  meta <- obj$meta
  time_chr <- as.character(meta[[p$time_col]])
  time_levels <- p$time_levels
  if (is.null(time_levels))
    time_levels <- .moduleviz_time_levels(meta[[p$time_col]])
  if (!is.null(p$group_col)) {
    group_chr <- as.character(meta[[p$group_col]])
    lookup <- unique(data.table(ckey = paste(group_chr, time_chr, sep = "|"),
                                Group = group_chr, Time = time_chr))
  } else {
    lookup <- unique(data.table(ckey = time_chr, Time = time_chr))
  }
  long[, key := as.character(key)]
  long <- merge(long, lookup, by.x = "key", by.y = "ckey", sort = FALSE)
  long[, Time := factor(Time, levels = time_levels)]

  if (!is.null(p$group_col)) {
    long[, Group := factor(Group, levels = unique(group_chr))]
    setorder(long, Group, Cluster, Time)
    group_cols <- setNames(.moduleviz_discrete_palette(nlevels(long$Group), palette),
                           levels(long$Group))
    ## explicit colours win for the levels they name
    if (!is.null(group_colors)) {
      named <- intersect(names(group_colors), names(group_cols))
      if (!length(named))
        warning("none of `group_colors` match the group levels: ",
                paste(levels(long$Group), collapse = ", "))
      group_cols[named] <- group_colors[named]
    }
    gg <- ggplot(long, aes(Time, z, group = Group, color = Group)) +
      scale_color_manual(values = group_cols)
  } else {
    setorder(long, Cluster, Time)
    cluster_cols <- setNames(
      .moduleviz_discrete_palette(nlevels(obj$cluster), cluster_palette),
      levels(obj$cluster))
    gg <- ggplot(long, aes(Time, z, group = Cluster, color = Cluster)) +
      scale_color_manual(values = cluster_cols, guide = "none")
  }

  ## One column of modules by default (strips on the right, like the heatmap's
  ## row titles); `ncol` > 1 wraps the modules into a grid.
  ncol <- max(1L, as.integer(ncol))
  n_facet_col <- min(ncol, nlevels(obj$cluster))
  gg <- gg +
    facet_wrap(~ ClusterF, ncol = n_facet_col,
               strip.position = if (n_facet_col == 1) "right" else "top") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_line(linewidth = line_width) + geom_point(size = point_size) +
    .moduleviz_theme(base_size) +
    theme(strip.text.y.right = element_text(angle = 0, hjust = 0),
          legend.title = element_blank(),
          plot.title = element_text(size = title_fontsize),
          axis.text.x = element_text(angle = x_angle,
                                     hjust = if (x_angle == 0) 0.5 else 1)) +
    labs(x = if (is.null(xlab)) p$time_col else xlab,
         y = "Mean z-score", color = p$group_col, title = title)

  ## optional relabelling of the ordering axis (e.g. 1/2/3 -> naive/mem/exh)
  if (!is.null(time_labels))
    gg <- gg + scale_x_discrete(labels = time_labels)

  if (!is.null(file)) {
    n_facet_row <- ceiling(nlevels(obj$cluster) / n_facet_col)
    if (is.null(width))
      width <- (if (!is.null(p$group_col)) 4.6 else 3.4) + 2.2 * (n_facet_col - 1)
    if (is.null(height))
      height <- (if (n_facet_col == 1) 1.0 else 1.5) * n_facet_row + 0.7
    ggsave(file, gg, width = width, height = height, limitsize = FALSE)
    message("wrote ", file)
  }
  gg
}

## =============================================================================
## Membership and summary export
## =============================================================================

#' Write the cluster memberships (and, optionally, centroids + stability).
#'
#' @param obj a "longi" or "multi_longi" object.
#' @param file path for the membership table (tab-separated).
#' @param prefix if given (and file is NULL), files are written as
#'   <prefix>_memberships.txt / _centroids.txt / _stability.txt.
#' @param write_extras also write centroids and stability tables.
#' @return invisibly, the membership data.table.
write_memberships <- function(obj, file = NULL, prefix = NULL,
                              write_extras = TRUE) {
  stopifnot(inherits(obj, c("longi", "multi_longi")))
  if (is.null(file)) {
    if (is.null(prefix)) stop("Supply either `file` or `prefix`.")
    file <- paste0(prefix, "_memberships.txt")
  }
  fwrite(obj$membership, file, sep = "\t", quote = FALSE)
  message("wrote ", file)

  if (write_extras) {
    base <- if (!is.null(prefix)) prefix else sub("\\.[^.]*$", "", file)
    if (!is.null(obj$centroids))
      fwrite(obj$centroids, paste0(base, "_centroids.txt"), sep = "\t", quote = FALSE)
    if (!is.null(obj$stability))
      fwrite(obj$stability, paste0(base, "_stability.txt"), sep = "\t", quote = FALSE)
  }
  invisible(obj$membership)
}

## =============================================================================
## Multi-omics heatmap
## =============================================================================

#' Collapse a feature x sample matrix to one row per gene.
#'
#'   - aggregate = "mean": average across that gene's features (NAs ignored);
#'   - aggregate = "top" : keep the single feature with the strongest |z|.
#'
#' @param z numeric matrix, rows = features.
#' @param genes gene key of each row.
#' @param aggregate "mean" or "top".
#' @return matrix with one row per gene, rownames = gene keys.
#' @noRd
.moduleviz_collapse_by_gene <- function(z, genes, aggregate = "mean") {
  keep  <- !is.na(genes) & nzchar(genes)
  z     <- z[keep, , drop = FALSE]
  genes <- genes[keep]
  if (aggregate == "mean") {
    ok  <- !is.na(z)
    out <- rowsum(ifelse(ok, z, 0), genes, reorder = FALSE) /
           rowsum(ok * 1, genes, reorder = FALSE)
  } else {
    score <- suppressWarnings(apply(abs(z), 1, max, na.rm = TRUE))
    o     <- order(genes, -score)
    o     <- o[!duplicated(genes[o])]
    out   <- z[o, , drop = FALSE]
    rownames(out) <- genes[o]
  }
  out
}

#' Normalise one entry of `layers` into a common internal description.
#'
#' A "longi" layer keeps its own z-scores, metadata, and column settings.  A raw
#' matrix is matched to its metadata (if any), its columns are ordered by group
#' then time, and it is z-scored later, once its final columns are known.
#' @noRd
.moduleviz_layer_spec <- function(x, name, meta, sample_col, time_col,
                                  group_col) {
  if (inherits(x, "longi")) {
    p <- x$params
    return(list(obj = x, mat = NULL, z = x$z, meta = x$meta,
                sample_col = p$sample_col, time_col = p$time_col,
                group_col = p$group_col))
  }
  ## a SummarizedExperiment layer brings its own sample metadata
  if (methods::is(x, "SummarizedExperiment")) {
    layer_meta <- if (is.null(meta) || is.data.frame(meta)) meta
                  else meta[[name]]
    input <- .moduleviz_unpack(x, layer_meta, sample_col = sample_col)
    x <- input$mat
    meta <- input$meta
    sample_col <- input$sample_col
  }
  m <- as.matrix(x)
  if (is.null(rownames(m)))
    stop(sprintf("layer '%s' has no rownames; supply feature IDs as rownames.",
                 name), call. = FALSE)
  lm <- if (is.null(meta) || is.data.frame(meta)) meta else meta[[name]]
  spec <- list(obj = NULL, mat = m, z = NULL, meta = NULL, raw_meta = lm,
               sample_col = sample_col, time_col = time_col,
               group_col = group_col)
  if (!is.null(lm)) {
    lm <- as.data.table(lm)
    if (!sample_col %in% names(lm))
      stop(sprintf("metadata for layer '%s' has no '%s' column.", name,
                   sample_col), call. = FALSE)
    lm <- lm[as.character(lm[[sample_col]]) %in% colnames(m)]
    if (nrow(lm)) {
      ord_cols <- intersect(c(group_col, time_col), names(lm))
      if (length(ord_cols))
        lm <- lm[do.call(order, c(as.list(lm[, ord_cols, with = FALSE]),
                                  list(seq_len(nrow(lm)))))]
      spec$mat  <- m[, as.character(lm[[sample_col]]), drop = FALSE]
      spec$meta <- lm
    }
  }
  if (!is.null(group_col) && !group_col %in% names(spec$meta))
    spec$group_col <- NULL
  spec
}

#' Cluster several omics layers into shared temporal modules, linked by gene.
#'
#' Generalises the paired RNA + ATAC view to any number of layers, e.g.
#' RNA + ATAC + H3K27ac + H3K27me3.  Every layer is z-scored within itself and
#' collapsed to one row per gene, so layers with very different feature types
#' (genes, peaks, CUT&RUN regions) line up row for row.
#'
#' Two clustering strategies:
#'   - `cluster_by = <layer>` (default: the first layer): that layer alone
#'     defines the modules and row order, exactly like `dual_heatmap()`;
#'     the other layers are shown in the same order.
#'   - `cluster_by = "joint"`: the per-gene z-scores of all layers are
#'     concatenated and clustered together, so a module groups genes whose
#'     expression AND chromatin trajectories agree.
#'
#' @param layers named list of layers.  Each entry is a "longi" object or a raw
#'   features x samples matrix.  Names become panel titles.
#' @param meta sample metadata for raw-matrix layers: one data.frame shared by
#'   all layers, or a named list with one data.frame per layer.  Ignored for
#'   "longi" layers, which carry their own.
#' @param sample_col,time_col,group_col metadata columns for raw-matrix layers
#'   (see longitudinal_cluster()).
#' @param cluster_by name or index of the layer that defines the modules, or
#'   "joint" to cluster all layers together.
#' @param gene_of function mapping row IDs to gene symbols, either one function
#'   for every layer or a list with one per layer (by name or position).  The
#'   default uses the first layer's IDs as-is and takes the text after the last
#'   "_" for every other layer ("PeakID_SYMBOL" -> "SYMBOL").
#' @param aggregate how to collapse several features per gene ("mean" or "top").
#' @param gene_set "intersect" keeps genes present in every layer; "reference"
#'   keeps every row of the `cluster_by` layer and leaves missing layers NA.
#' @param balance_layers joint clustering only: weight each layer's columns by
#'   1/sqrt(number of columns) so a layer with many samples does not dominate.
#' @param ... passed to longitudinal_cluster() (k, method, stability, ...) when
#'   the clustering layer is a raw matrix, or always for "joint".
#' @return object of class "multi_longi": a list with layers (per-layer aligned
#'   z matrix and metadata), cluster, genes, membership, stability, reference
#'   (the "longi" object the modules came from), params.
multi_cluster <- function(layers, meta = NULL,
                                sample_col = "SampleID",
                                time_col   = "TimeHr",
                                group_col  = NULL,
                                cluster_by = 1,
                                gene_of = NULL,
                                aggregate = c("mean", "top"),
                                gene_set = c("intersect", "reference"),
                                balance_layers = TRUE,
                                ...) {
  aggregate <- match.arg(aggregate)
  gene_set  <- match.arg(gene_set)
  if (inherits(layers, "longi") || !is.list(layers) || length(layers) < 2)
    stop("`layers` must be a list of at least two layers.", call. = FALSE)
  nms <- names(layers)
  if (is.null(nms)) nms <- rep("", length(layers))
  nms[!nzchar(nms)] <- paste0("Layer", which(!nzchar(nms)))
  nms <- make.unique(nms)
  names(layers) <- nms
  n_layer <- length(layers)

  ## ---- which layer drives the modules -------------------------------------
  joint <- identical(cluster_by, "joint")
  if (!joint) {
    ref_i <- if (is.character(cluster_by)) match(cluster_by, nms) else as.integer(cluster_by)
    if (length(ref_i) != 1 || is.na(ref_i) || ref_i < 1 || ref_i > n_layer)
      stop("`cluster_by` must be a layer name, a layer index, or \"joint\".",
           call. = FALSE)
  } else if (gene_set == "reference") {
    stop("gene_set = \"reference\" needs a reference layer; joint clustering ",
         "uses the genes shared by every layer.", call. = FALSE)
  }

  ## ---- gene mappers: one per layer ----------------------------------------
  default_gene_of <- function(i) if (i == 1) function(id) id
                                 else function(id) sub("^.*_", "", id)
  if (is.function(gene_of)) gene_of <- rep(list(gene_of), n_layer)
  if (is.null(gene_of)) gene_of <- list()
  gene_of <- lapply(seq_len(n_layer), function(i) {
    f <- if (!is.null(names(gene_of)) && nms[i] %in% names(gene_of)) gene_of[[nms[i]]]
         else if (is.null(names(gene_of)) && i <= length(gene_of)) gene_of[[i]]
         else NULL
    if (is.null(f)) default_gene_of(i) else f
  })

  specs <- Map(.moduleviz_layer_spec, layers, nms,
               MoreArgs = list(meta = meta, sample_col = sample_col,
                               time_col = time_col, group_col = group_col))

  if (!joint) {
    ## ---- reference mode: one layer clusters, the rest follow ---------------
    ref <- specs[[ref_i]]
    if (is.null(ref$obj)) {
      if (is.null(ref$raw_meta))
        stop(sprintf("layer '%s' is a raw matrix; supply `meta` so it can be clustered.",
                     nms[ref_i]), call. = FALSE)
      ref$obj <- longitudinal_cluster(ref$mat, ref$raw_meta,
                                      sample_col = ref$sample_col,
                                      time_col = time_col,
                                      group_col = group_col, ...)
      specs[[ref_i]] <- .moduleviz_layer_spec(ref$obj, nms[ref_i])
    } else if (length(list(...))) {
      warning("clustering arguments in `...` are ignored: layer '", nms[ref_i],
              "' is already clustered.", call. = FALSE)
    }
    ref_obj <- ref$obj
    ref_z   <- ref_obj$z
    g_ref   <- gene_of[[ref_i]](rownames(ref_z))
  }

  ## Raw layers without metadata follow the clustering layer's sample order
  ## where the column names overlap.
  anchor_cols <- if (!joint) colnames(ref_z) else NULL
  for (i in seq_len(n_layer)) {
    s <- specs[[i]]
    if (!is.null(s$z)) next
    if (is.null(s$meta) && !is.null(anchor_cols)) {
      common <- intersect(anchor_cols, colnames(s$mat))
      if (length(common) >= 2) s$mat <- s$mat[, common, drop = FALSE]
    }
    s$z <- zscore_rows(s$mat)
    specs[[i]] <- s
  }

  collapse_layer <- function(i)
    .moduleviz_collapse_by_gene(specs[[i]]$z,
                                gene_of[[i]](rownames(specs[[i]]$z)),
                                aggregate)

  if (!joint) {
    others    <- setdiff(seq_len(n_layer), ref_i)
    collapsed <- setNames(vector("list", n_layer), nms)
    for (i in others) collapsed[[i]] <- collapse_layer(i)
    in_layer  <- vapply(others, function(i) g_ref %in% rownames(collapsed[[i]]),
                        logical(length(g_ref)))
    in_layer  <- matrix(in_layer, nrow = length(g_ref))
    keep_rows <- if (gene_set == "intersect") which(rowSums(in_layer) == length(others))
                 else seq_along(g_ref)
    if (!length(keep_rows))
      stop("No genes shared between the layers; check the `gene_of` mappers.",
           call. = FALSE)

    row_ids <- rownames(ref_z)[keep_rows]
    genes   <- g_ref[keep_rows]
    cluster <- droplevels(ref_obj$cluster[keep_rows])
    layer_z <- setNames(vector("list", n_layer), nms)
    layer_z[[ref_i]] <- ref_z[keep_rows, , drop = FALSE]
    for (i in others) {
      zi <- collapsed[[i]][match(genes, rownames(collapsed[[i]])), , drop = FALSE]
      rownames(zi) <- row_ids
      layer_z[[i]] <- zi
    }
  } else {
    ## ---- joint mode: cluster the concatenated per-gene profiles ------------
    if (isTRUE(list(...)$aggregate_replicates))
      stop("aggregate_replicates is not supported for joint clustering; pass ",
           "per-layer longitudinal_cluster(..., aggregate_replicates = TRUE) ",
           "objects instead.", call. = FALSE)
    no_time <- nms[vapply(specs, function(s)
      is.null(s$meta) || !s$time_col %in% names(s$meta), logical(1))]
    if (length(no_time))
      stop("joint clustering needs sample metadata with a time column for every ",
           "layer; missing for: ", paste(no_time, collapse = ", "), call. = FALSE)

    collapsed <- setNames(lapply(seq_len(n_layer), collapse_layer), nms)
    genes <- Reduce(intersect, lapply(collapsed, rownames))
    if (length(genes) < 2)
      stop("Fewer than 2 genes shared by every layer; check the `gene_of` mappers.",
           call. = FALSE)

    ## Shared time / group column names come from the first layer; the group
    ## is used only when every layer has one.
    tc <- specs[[1]]$time_col
    gc <- if (all(vapply(specs, function(s) !is.null(s$group_col), logical(1))))
            specs[[1]]$group_col else NULL
    blocks <- lapply(seq_len(n_layer), function(i) {
      zi <- collapsed[[i]][genes, , drop = FALSE]
      w  <- if (balance_layers) 1 / sqrt(ncol(zi)) else 1
      colnames(zi) <- paste0(nms[i], "::", colnames(zi))
      zi * w
    })
    joint_meta <- rbindlist(lapply(seq_len(n_layer), function(i) {
      s  <- specs[[i]]
      m  <- s$meta[match(colnames(s$z), as.character(s$meta[[s$sample_col]]))]
      dt <- data.table(LayerSample = paste0(nms[i], "::", colnames(s$z)),
                       Layer = nms[i])
      dt[[tc]] <- m[[s$time_col]]
      if (!is.null(gc)) dt[[gc]] <- as.character(m[[s$group_col]])
      dt
    }))
    ref_obj <- longitudinal_cluster(do.call(cbind, blocks), joint_meta,
                                    sample_col = "LayerSample", time_col = tc,
                                    group_col = gc, is_zscore = TRUE, ...)

    ## Order rows within each module by the unweighted signal.
    genes   <- rownames(ref_obj$z)
    layer_z <- lapply(collapsed, function(zi) zi[genes, , drop = FALSE])
    score   <- suppressWarnings(apply(abs(do.call(cbind, layer_z)), 1, max,
                                      na.rm = TRUE))
    o       <- order(ref_obj$cluster, -score)
    genes   <- genes[o]
    layer_z <- lapply(layer_z, function(zi) zi[o, , drop = FALSE])
    cluster <- ref_obj$cluster[o]
    row_ids <- genes
    ref_i   <- NA_integer_
  }

  out_layers <- setNames(lapply(seq_len(n_layer), function(i) {
    s <- specs[[i]]
    list(z = layer_z[[i]], meta = s$meta, sample_col = s$sample_col,
         time_col = s$time_col, group_col = s$group_col)
  }), nms)

  score <- suppressWarnings(apply(abs(do.call(cbind, layer_z)), 1, max,
                                  na.rm = TRUE))
  membership <- data.table(ID = row_ids, Gene = genes,
                           Cluster = as.character(cluster),
                           MaxAbsZ = round(score, 3))

  structure(list(
    layers     = out_layers,
    cluster    = cluster,
    genes      = genes,
    membership = membership,
    stability  = ref_obj$stability,
    reference  = ref_obj,
    params     = list(cluster_by = if (joint) "joint" else nms[ref_i],
                      layer_names = nms, aggregate = aggregate,
                      gene_set = gene_set, balance_layers = balance_layers)
  ), class = "multi_longi")
}

print.multi_longi <- function(x, ...) {
  p <- x$params
  cat(sprintf("<multi_longi>  %d genes x %d layers (%s),  %d modules\n",
              length(x$genes), length(x$layers),
              paste(p$layer_names, collapse = " + "), nlevels(x$cluster)))
  cat(sprintf("  clustered by: %s\n", p$cluster_by))
  if (!is.null(x$stability))
    cat(sprintf("  stability (mean Jaccard): %.2f - %.2f\n",
                min(x$stability$MeanJaccard), max(x$stability$MeanJaccard)))
  invisible(x)
}

#' Draw several omics layers side by side, linked by gene.
#'
#' One panel per layer (e.g. RNA | ATAC | H3K27ac | H3K27me3), all sharing the
#' module split and row order from multi_cluster().  The layers share one
#' z-score colour scale, so a single legend covers every panel.  Each panel gets
#' its own metadata annotation bars, coloured consistently across layers.
#'
#' @param x a "multi_longi" object, or a named list of layers; a list is passed
#'   to multi_cluster() together with `...`.
#' @param annotation_cols metadata columns to show as top bars on each layer.
#'   Defaults to each layer's group_col (plus time_col when the time points
#'   are not printed under the columns); `character(0)` for none.
#' @param annotation_colors optional named list mapping annotation -> colours.
#' @param label_features genes / IDs to label on the right of the last panel.
#' @param top_n_label per-module automatic labels, ranked by the largest |z|
#'   across all layers (0 = none).
#' @param col_fun z-score colour map.  Overrides `zscore_palette` when supplied.
#' @param zscore_palette "rdbu" (default) or "viridis", shared by every layer.
#' @param zscore_limits numeric length-2 z-score saturation range.
#' @param cluster_palette discrete palette for the module annotation bar.
#' @param annotation_palette discrete palette for categorical metadata bars.
#' @param na_col colour for genes missing from a layer (gene_set = "reference").
#' @param label_fontsize,row_title_fontsize,column_name_fontsize text sizes.
#' @param annotation_name_fontsize,legend_fontsize,title_fontsize text sizes for
#'   the annotation names, legends, and panel / figure titles.
#' @param show_column_names show sample names along the bottom of every layer.
#' @param layer_gap gap between the layer panels, in mm.
#' @param column_labels "time" (default) prints only the time point under each
#'   column; "sample" prints the full sample names.
#' @param layer_widths "equal" (default) gives every layer the same panel
#'   width; "samples" makes each panel as wide as its number of samples.
#' @param title optional overall title; none by default (each panel is already
#'   headed by its layer name).
#' @param file optional PDF or PNG path.
#' @param width,height size in inches; width defaults to 2.5 + 2 per layer.
#' @param ... passed to multi_cluster() when `x` is a list of layers.
#' @return invisibly, the "multi_longi" object with the drawn HeatmapList in
#'   `$heatmap`.
multi_heatmap <- function(x,
                                annotation_cols = NULL,
                                annotation_colors = NULL,
                                label_features = NULL,
                                top_n_label = 0,
                                col_fun = NULL,
                                zscore_palette = c("rdbu", "viridis"),
                                zscore_limits = c(-2, 2),
                                cluster_palette = "paired",
                                annotation_palette = "set2",
                                na_col = "grey85",
                                label_fontsize = 9,
                                row_title_fontsize = 10,
                                column_name_fontsize = 9,
                                annotation_name_fontsize = 9,
                                legend_fontsize = 9,
                                title_fontsize = 11,
                                show_column_names = TRUE,
                                layer_gap = 5,
                                title = NULL,
                                file = NULL, width = NULL, height = 6.5, ...,
                                enrichment = NULL,
                                enrichment_top_n = 3,
                                enrichment_padj = 0.05,
                                enrichment_fontsize = 9,
                                column_labels = c("time", "sample"),
                                layer_widths = c("equal", "samples")) {
  column_labels <- match.arg(column_labels)
  layer_widths <- match.arg(layer_widths)
  if (!inherits(x, "multi_longi")) x <- multi_cluster(x, ...)
  if (is.null(col_fun))
    col_fun <- .moduleviz_zscore_col_fun(zscore_palette, zscore_limits)
  layers  <- x$layers
  nms     <- names(layers)
  n_layer <- length(layers)
  cluster <- x$cluster
  if (is.null(width)) width <- 2.5 + 2 * n_layer
  if (isFALSE(title)) title <- NULL
  clabs <- lapply(layers, function(L)
    .moduleviz_column_labels(L$z, L$meta, L$time_col, column_labels))

  ## ---- annotation bars, coloured consistently across layers ---------------
  ann_by_layer <- lapply(layers, function(L) {
    cols <- if (is.null(annotation_cols)) {
      time_on_axis <- show_column_names && .moduleviz_column_labels(
        L$z, L$meta, L$time_col, column_labels)$is_time
      c(L$group_col, if (!time_on_axis) L$time_col)
    } else annotation_cols
    intersect(cols, names(L$meta))
  })
  all_ann <- unique(unlist(ann_by_layer))
  missing_anno <- setdiff(annotation_cols, all_ann)
  if (length(missing_anno))
    warning("annotation_cols not found in any layer's metadata and skipped: ",
            paste(missing_anno, collapse = ", "), call. = FALSE)
  col_list <- if (is.null(annotation_colors)) list() else annotation_colors
  for (a in all_ann) {
    if (!is.null(col_list[[a]])) next
    vals <- lapply(layers[vapply(ann_by_layer, function(s) a %in% s, logical(1))],
                   function(L) L$meta[[a]])
    if (all(vapply(vals, is.numeric, logical(1)))) {
      col_list[[a]] <- circlize::colorRamp2(range(unlist(vals), na.rm = TRUE),
                                            c("#deebf7", "#08519c"))
    } else {
      lv <- unique(unlist(lapply(vals, function(v)
        if (is.factor(v)) intersect(levels(v), as.character(v))
        else unique(as.character(v)))))
      col_list[[a]] <- setNames(
        .moduleviz_discrete_palette(length(lv), annotation_palette), lv)
    }
  }
  ## each legend once (first layer that has it), each name once (last layer)
  has_ann    <- lapply(all_ann, function(a)
    which(vapply(ann_by_layer, function(s) a %in% s, logical(1))))
  first_with <- setNames(vapply(has_ann, min, numeric(1)), all_ann)
  last_with  <- setNames(vapply(has_ann, max, numeric(1)), all_ann)

  ## ---- row titles, module bar, labels --------------------------------------
  cn <- table(cluster)
  row_title <- sprintf("%s\n(n=%d)", levels(cluster), as.integer(cn))
  cluster_cols <- .moduleviz_cluster_cols(levels(cluster), cluster_palette)
  left_anno <- rowAnnotation(
    Cluster = cluster, col = list(Cluster = cluster_cols),
    show_annotation_name = FALSE, show_legend = FALSE,
    simple_anno_size = unit(2.5, "mm"))

  at <- integer(0)
  if (length(label_features))
    at <- which(x$genes %in% label_features | x$membership$ID %in% label_features)
  if (top_n_label > 0) {
    dt <- data.table(idx = seq_along(x$genes), Cluster = x$membership$Cluster,
                     value = x$membership$MaxAbsZ)
    setorder(dt, -value)
    at <- c(at, dt[, head(.SD, top_n_label), by = Cluster]$idx)
  }
  at <- sort(unique(at))
  right_parts <- list()
  if (length(at))
    right_parts$mark <- anno_mark(which = "row",
      at = at, labels = x$genes[at],
      labels_gp = gpar(fontsize = label_fontsize, fontface = "italic"),
      link_gp = gpar(lwd = 0.5), link_width = unit(4, "mm"))
  if (!is.null(enrichment))
    right_parts$Enrichment <- .moduleviz_enrichment_anno(
      enrichment, levels(cluster), top_n = enrichment_top_n,
      padj_cutoff = enrichment_padj, fontsize = enrichment_fontsize)
  right_anno <- NULL
  if (length(right_parts))
    right_anno <- do.call(rowAnnotation,
                          c(right_parts, list(show_annotation_name = FALSE)))

  ## ---- one panel per layer -------------------------------------------------
  panel <- function(i) {
    L <- layers[[i]]
    top_anno <- NULL
    a <- ann_by_layer[[i]]
    if (length(a)) {
      top_anno <- HeatmapAnnotation(
        df = as.data.frame(L$meta[, a, with = FALSE]), col = col_list[a],
        show_legend = unname(first_with[a] == i),
        show_annotation_name = unname(last_with[a] == i),
        simple_anno_size = unit(3, "mm"),
        annotation_name_gp = gpar(fontsize = annotation_name_fontsize),
        annotation_legend_param = .moduleviz_legend_gp(legend_fontsize))
    }
    col_split <- if (!is.null(L$group_col) && L$group_col %in% names(L$meta))
      factor(L$meta[[L$group_col]], levels = unique(L$meta[[L$group_col]]))
    first <- i == 1
    Heatmap(
      L$z, name = nms[i], col = col_fun, na_col = na_col,
      width = if (layer_widths == "equal") unit(1, "null") else unit(ncol(L$z), "null"),
      row_split = if (first) cluster, cluster_rows = FALSE,
      cluster_row_slices = FALSE,
      row_title = if (first) row_title, row_title_rot = 0,
      row_title_gp = gpar(fontsize = row_title_fontsize),
      row_gap = unit(1, "mm"), show_row_names = FALSE,
      column_split = col_split, cluster_columns = FALSE,
      column_title = nms[i], column_gap = unit(1, "mm"),
      column_title_gp = gpar(fontsize = title_fontsize),
      show_column_names = show_column_names,
      column_labels = clabs[[i]]$labels, column_names_rot = clabs[[i]]$rot,
      column_names_centered = clabs[[i]]$rot == 0,
      column_names_gp = gpar(fontsize = column_name_fontsize),
      column_names_max_height = max_text_width(
        clabs[[i]]$labels, gp = gpar(fontsize = column_name_fontsize)),
      top_annotation = top_anno,
      left_annotation = if (first) left_anno,
      right_annotation = if (i == n_layer) right_anno,
      ## the layers share one colour scale, so one legend covers them all
      show_heatmap_legend = first,
      use_raster = TRUE, raster_quality = 2,
      heatmap_legend_param = c(list(title = "z-score", direction = "horizontal",
                                    title_position = "topcenter",
                                    legend_width = unit(2.5, "cm"),
                                    grid_height = unit(2.5, "mm")),
                               .moduleviz_legend_gp(legend_fontsize)))
  }
  combined <- Reduce(`+`, lapply(seq_len(n_layer), panel))

  draw_it <- function() draw(combined, heatmap_legend_side = "bottom",
                             annotation_legend_side = "bottom",
                             merge_legend = TRUE, ht_gap = unit(layer_gap, "mm"),
                             column_title = title,
                             column_title_gp = gpar(fontsize = title_fontsize))
  if (!is.null(file)) {
    .moduleviz_write_plot(file, width, height, draw_it)
  } else {
    draw_it()
  }
  x$heatmap <- combined
  invisible(x)
}

#' Module trajectories across layers: one line per layer, one facet per module.
#'
#' For each layer, every sample contributes its mean z-score over the module's
#' genes; points are the mean of those values at each time point and the error
#' bars their SD across replicates.  Answers "does the chromatin mark rise and
#' fall with expression in this module?" at a glance.
#'
#' @param x a "multi_longi" object (or the value of multi_heatmap()).
#' @param file optional output path (any ggsave() format).
#' @param width,height size in inches; auto-scaled when NULL.
#' @param ncol facet columns when the layers have no group_col.
#' @param palette discrete palette for the layer lines (default "okabe_ito").
#' @param layer_colors optional named vector of explicit layer colours;
#'   overrides `palette` for the layers it names.
#' @param error_bars draw +/- SD across replicates where there are replicates.
#' @param base_size,title,title_fontsize,xlab,time_labels,x_angle,line_width,point_size
#'   as in pattern_lineplot().
#' @return a ggplot object (also saved if `file` is given).
multi_lineplot <- function(x, file = NULL, width = NULL, height = NULL,
                                 ncol = 1,
                                 split_groups = FALSE,
                                 palette = "okabe_ito",
                                 layer_colors = NULL,
                                 error_bars = TRUE,
                                 base_size = 10,
                                 title = NULL,
                                 title_fontsize = NULL,
                                 xlab = NULL,
                                 time_labels = NULL,
                                 x_angle = 0,
                                 line_width = 0.7,
                                 point_size = 1.5) {
  stopifnot(inherits(x, "multi_longi"))
  if (is.null(title_fontsize)) title_fontsize <- base_size
  cl_levels <- levels(x$cluster)

  ## per-sample module means, one row per layer x module x sample
  has_time <- vapply(x$layers, function(L)
    !is.null(L$meta) && L$time_col %in% names(L$meta), logical(1))
  if (any(!has_time))
    warning("layers without a time column in their metadata are left out: ",
            paste(names(x$layers)[!has_time], collapse = ", "), call. = FALSE)
  if (!any(has_time)) stop("No layer has sample metadata with a time column.")
  used <- names(x$layers)[has_time]

  ## Time points across layers: numeric when every layer is numeric, otherwise
  ## each layer's own ordering, merged in layer order.
  time_vals <- lapply(x$layers[used], function(L) L$meta[[L$time_col]])
  time_levels <- if (all(vapply(time_vals, is.numeric, logical(1)))) {
    as.character(sort(unique(unlist(time_vals))))
  } else {
    unique(unlist(lapply(time_vals, .moduleviz_time_levels)))
  }

  per <- rbindlist(lapply(used, function(nm) {
    L <- x$layers[[nm]]
    m   <- L$meta[match(colnames(L$z), as.character(L$meta[[L$sample_col]]))]
    grp <- if (!is.null(L$group_col) && L$group_col %in% names(m))
             as.character(m[[L$group_col]]) else NA_character_
    rbindlist(lapply(cl_levels, function(cl)
      data.table(Layer = nm, Cluster = cl, Group = grp,
                 Time = as.character(m[[L$time_col]]),
                 z = colMeans(L$z[x$cluster == cl, , drop = FALSE], na.rm = TRUE))))
  }))

  summ <- per[, .(z = mean(z, na.rm = TRUE), zsd = sd(z, na.rm = TRUE)),
              by = .(Layer, Cluster, Group, Time)]
  summ[, Time := factor(Time, levels = time_levels)]
  summ[, Layer := factor(Layer, levels = intersect(names(x$layers), Layer))]
  cn <- table(x$cluster)
  facet_lab <- setNames(sprintf("%s (n=%d)", cl_levels, as.integer(cn)), cl_levels)
  summ[, ClusterF := factor(facet_lab[Cluster], levels = facet_lab)]

  layer_cols <- setNames(.moduleviz_discrete_palette(nlevels(summ$Layer), palette),
                         levels(summ$Layer))
  if (!is.null(layer_colors)) {
    named <- intersect(names(layer_colors), names(layer_cols))
    layer_cols[named] <- layer_colors[named]
  }

  grouped <- any(!is.na(summ$Group))
  ## Default layout: one column of modules.  Groups are then told apart by
  ## line type; `split_groups = TRUE` gives every group its own column instead.
  by_linetype <- grouped && !split_groups
  if (grouped) {
    summ[, Group := factor(Group, levels = unique(Group))]
    summ[, LineGroup := interaction(Layer, Group, drop = TRUE)]
  } else {
    summ[, LineGroup := Layer]
  }

  pd <- position_dodge(width = 0.15)
  gg <- ggplot(summ, aes(Time, z, group = LineGroup, color = Layer)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75")
  gg <- gg + if (by_linetype) {
    geom_line(aes(linetype = Group), linewidth = line_width, position = pd)
  } else {
    geom_line(linewidth = line_width, position = pd)
  }
  if (error_bars && any(!is.na(summ$zsd)))
    gg <- gg + geom_errorbar(data = summ[!is.na(zsd)],
                             aes(ymin = z - zsd, ymax = z + zsd),
                             width = 0, linewidth = 0.4, position = pd)
  gg <- gg + geom_point(size = point_size, position = pd) +
    scale_color_manual(values = layer_cols)

  if (grouped && split_groups) {
    gg <- gg + facet_grid(rows = vars(ClusterF), cols = vars(Group))
    n_col <- nlevels(summ$Group)
    n_row <- length(cl_levels)
  } else {
    n_col <- min(max(1L, as.integer(ncol)), length(cl_levels))
    gg <- gg + facet_wrap(~ ClusterF, ncol = n_col,
                          strip.position = if (n_col == 1) "right" else "top")
    n_row <- ceiling(length(cl_levels) / n_col)
  }

  gg <- gg +
    .moduleviz_theme(base_size) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          strip.text.y.right = element_text(angle = 0, hjust = 0),
          plot.title = element_text(size = title_fontsize),
          axis.text.x = element_text(angle = x_angle,
                                     hjust = if (x_angle == 0) 0.5 else 1)) +
    labs(x = if (is.null(xlab)) x$layers[[1]]$time_col else xlab,
         y = "Mean z-score", title = title)
  if (!is.null(time_labels))
    gg <- gg + scale_x_discrete(labels = time_labels)

  if (!is.null(file)) {
    if (is.null(width))  width  <- 2.2 * n_col + 1.6
    if (is.null(height)) height <- (if (n_col == 1) 1.2 else 1.6) * n_row + 1.0
    ggsave(file, gg, width = width, height = height, limitsize = FALSE)
    message("wrote ", file)
  }
  gg
}

## =============================================================================
## Paired omics heatmap
## =============================================================================

#' Draw two omics layers side by side, linked by gene.
#'
#' The two-layer case of multi_heatmap(), kept for its simpler arguments:
#' a primary layer (e.g. RNA, one row per gene) next to a secondary layer (e.g.
#' ATAC / CUT&RUN peaks), sharing a common set of genes and row order.
#'
#' Because one gene can map to several peaks, the secondary layer is collapsed to
#' one value per gene:
#'   - aggregate = "mean": average z-score across that gene's peaks;
#'   - aggregate = "top" : the single peak with the strongest |z|.
#'
#' Genes are clustered/ordered using the *primary* layer (so RNA drives the
#' module structure), then the secondary layer is shown in the same order.
#'
#' @param primary,secondary "longi" objects (from longitudinal_cluster) OR raw
#'   matrices with gene-bearing rownames.  If raw matrices are passed, provide
#'   `meta` and the *_col arguments; the primary is clustered here.
#' @param primary_name,secondary_name legend / title labels.
#' @param gene_of_primary,gene_of_secondary functions mapping a feature ID to a
#'   gene symbol.  Default assumes "PeakID_SYMBOL" (takes text after last "_")
#'   for the secondary and uses the ID as-is for the primary.
#' @param aggregate how to collapse secondary peaks per gene ("mean" or "top").
#' @param label_features genes to label on the right (character vector).
#' @param annotation_cols,annotation_colors top annotation bars and their
#'   colours, as in multi_heatmap(); by default the group column of each layer
#'   that carries sample metadata.
#' @param col_fun z-score colour map.  Overrides `zscore_palette` when supplied.
#' @param zscore_palette continuous palette for both layers: "rdbu" (default) or
#'   "viridis".  Both layers share one scale so they stay comparable.
#' @param zscore_limits numeric length-2 z-score saturation range.
#' @param cluster_palette discrete palette for the module annotation bar.
#'   Defaults to "paired".
#' @param label_fontsize,row_title_fontsize,column_name_fontsize text sizes.
#' @param legend_fontsize,title_fontsize text sizes for the legends and the
#'   per-layer column titles.
#' @param show_column_names show sample names along the bottom of both layers.
#' @param file optional PDF or PNG path (this is a *separate* new figure).
#' @param width,height PDF size.
#' @param ... passed through to longitudinal_cluster when raw matrices are given.
#' @return invisibly, a list(primary_obj, secondary_matrix, genes) plus the drawn
#'   combined heatmap.
dual_heatmap <- function(primary, secondary,
                               primary_name = "RNA", secondary_name = "ATAC",
                               gene_of_primary   = function(id) id,
                               gene_of_secondary = function(id) sub("^.*_", "", id),
                               aggregate = c("mean", "top"),
                               label_features = NULL,
                               col_fun = NULL,
                               zscore_palette = c("rdbu", "viridis"),
                               zscore_limits = c(-2, 2),
                               cluster_palette = "paired",
                               label_fontsize = 9,
                               row_title_fontsize = 10,
                               column_name_fontsize = 9,
                               legend_fontsize = 9,
                               title_fontsize = 11,
                               show_column_names = TRUE,
                               file = NULL, width = 6.5, height = 6.5, ...,
                               enrichment = NULL, enrichment_top_n = 3,
                               enrichment_padj = 0.05, annotation_cols = NULL,
                               annotation_colors = NULL) {
  aggregate <- match.arg(aggregate)

  ## Allow raw matrices: cluster the primary here.
  if (!inherits(primary, "longi")) {
    primary <- longitudinal_cluster(primary, ...)
  }

  x <- multi_cluster(
    setNames(list(primary, secondary), c(primary_name, secondary_name)),
    gene_of = list(gene_of_primary, gene_of_secondary), aggregate = aggregate)
  x <- multi_heatmap(
    x, annotation_cols = annotation_cols,
    annotation_colors = annotation_colors, label_features = label_features,
    col_fun = col_fun, zscore_palette = zscore_palette,
    zscore_limits = zscore_limits, cluster_palette = cluster_palette,
    label_fontsize = label_fontsize, row_title_fontsize = row_title_fontsize,
    column_name_fontsize = column_name_fontsize,
    legend_fontsize = legend_fontsize, title_fontsize = title_fontsize,
    show_column_names = show_column_names,
    file = file, width = width, height = height,
    enrichment = enrichment, enrichment_top_n = enrichment_top_n,
    enrichment_padj = enrichment_padj)

  invisible(list(primary = primary,
                 z_secondary = x$layers[[2]]$z,
                 genes = x$genes,
                 heatmap = x$heatmap))
}

## =============================================================================
## Convenience wrapper: run the whole single-omics pipeline in one call
## =============================================================================

#' One-call pipeline: cluster + heatmap + line plot + write memberships.
#'
#' @param mat,meta,... passed to longitudinal_cluster().
#' @param out_prefix path prefix for all output files.
#' @param label_features,top_n_label,annotation_cols passed to longitudinal_heatmap().
#' @param zscore_palette,cluster_palette,annotation_palette palettes forwarded to
#'   longitudinal_heatmap() / pattern_lineplot().
#' @param heatmap_width,heatmap_height heatmap PDF size.
#' @return the "longi" object (invisibly).
run_longitudinal <- function(mat, meta = NULL, out_prefix,
                             label_features = NULL, top_n_label = 5,
                             annotation_cols = NULL, rank_table = NULL,
                             zscore_palette = "rdbu",
                             cluster_palette = "paired",
                             annotation_palette = "set2",
                             heatmap_width = 5, heatmap_height = 6.5, ...) {
  obj <- longitudinal_cluster(mat, meta, ...)
  longitudinal_heatmap(obj, annotation_cols = annotation_cols,
                       label_features = label_features, top_n_label = top_n_label,
                       rank_table = rank_table,
                       zscore_palette = zscore_palette,
                       cluster_palette = cluster_palette,
                       annotation_palette = annotation_palette,
                       file = paste0(out_prefix, "_heatmap.pdf"),
                       width = heatmap_width, height = heatmap_height)
  pattern_lineplot(obj, file = paste0(out_prefix, "_pattern_lines.pdf"),
                   cluster_palette = cluster_palette)
  write_memberships(obj, prefix = out_prefix)
  invisible(obj)
}
