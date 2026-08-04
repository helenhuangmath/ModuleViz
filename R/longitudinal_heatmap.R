## =============================================================================
## longitudinal_heatmap.R
##
## A small toolkit of functions for longitudinal / time-course omics analysis
## (RNA-seq, ATAC-seq, Cut&Run, drug-treatment time courses, ...).
##
## Given a matrix of features (genes or peaks) measured across an ordered set of
## samples (time points / conditions), these functions:
##
##   1. cluster the features into temporal patterns / modules, with a choice of
##      algorithm (k-means, hierarchical, PAM) and a stability assessment so the
##      modules are reproducible;                              -> longitudinal_cluster()
##   2. draw an ordered heatmap running from "high at the first time point" to
##      "high at the last time point", with metadata colour bars;   -> longitudinal_heatmap()
##   3. draw a line/pattern plot of the mean trajectory per module;  -> pattern_lineplot()
##   4. write out the cluster memberships (and centroids, stability);-> write_memberships()
##   5. let the user pick specific features to label on the heatmap; -> label_features= arg
##   6. put two data sets side by side (e.g. RNA + ATAC, RNA + Cut&Run) so genes
##      and peaks are shown together in one figure.                  -> dual_omics_heatmap()
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

.moduleviz_write_plot <- function(file, width, height, draw_it) {
  ext <- tolower(sub("^.*\\.([^.]+)$", "\\1", file))
  if (identical(ext, "png")) {
    grDevices::png(file, width = width, height = height, units = "in", res = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_it()
  } else {
    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_it()
  }
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
  out <- sapply(group_levels, function(g) {
    cols <- which(grp == g)
    if (length(cols) == 0) return(rep(NA_real_, nrow(mat)))
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  out <- matrix(out, nrow = nrow(mat),
                dimnames = list(rownames(mat), group_levels))
  out
}

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
  method <- match.arg(method, c("kmeans", "hierarchical", "pam"))
  n_uniq <- nrow(unique(z))
  k <- min(k, n_uniq)
  if (k < 1) stop("Need at least one distinct row to cluster.")

  ## "correlation" distance = 1 - Pearson correlation between feature profiles.
  make_dist <- function(mat) {
    if (dist_method == "correlation") {
      as.dist(1 - cor(t(mat)))
    } else {
      dist(mat, method = dist_method)
    }
  }

  if (method == "kmeans") {
    set.seed(seed)
    ## Hartigan-Wong (R's default) rarely leaves empty clusters, which makes it
    ## more robust than Lloyd on tightly structured time-course data.
    km <- kmeans(z, centers = k, nstart = nstart, iter.max = 1000)
    return(km$cluster)
  }
  if (method == "hierarchical") {
    hc <- hclust(make_dist(z), method = hclust_method)
    return(cutree(hc, k = k))
  }
  ## pam (partitioning around medoids) -- robust k-means alternative
  if (!requireNamespace("cluster", quietly = TRUE))
    stop("method='pam' needs the 'cluster' package.")
  set.seed(seed)
  pm <- cluster::pam(make_dist(z), k = k, diss = TRUE)
  as.integer(pm$clustering)
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
  centers  <- t(sapply(labs, function(cl)
    colMeans(tp_mean[cluster == cl, , drop = FALSE], na.rm = TRUE)))
  colnames(centers) <- time_levels
  rownames(centers) <- labs

  peak_at <- apply(centers, 1, which.max)
  slope   <- centers[, ncol(centers)] - centers[, 1]
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
  n     <- nrow(z)
  labs  <- sort(unique(cluster))
  jacc  <- matrix(NA_real_, nrow = length(labs), ncol = n_resample,
                  dimnames = list(as.character(labs), NULL))

  for (b in seq_len(n_resample)) {
    set.seed(seed + b)
    idx <- sample.int(n, size = max(2L, floor(subsample_frac * n)))
    zb  <- z[idx, , drop = FALSE]
    ## reclustering the subsample; use a per-iteration seed so starts differ.
    ## Empty-cluster warnings are expected on small resamples and handled by
    ## kmeans internally, so silence them here to keep the console readable.
    cb  <- tryCatch(
      suppressWarnings(.cluster_once(zb, k = k, method = method, seed = seed + b,
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
    Cluster     = paste0("Cluster", labs),
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
                     chosen_k = NULL, annotate_k = TRUE, base_size = 11) {
  k_range <- k_range[k_range < nrow(unique(z))]
  d_for_sil <- if (dist_method == "correlation") as.dist(1 - cor(t(z)))
               else dist(z, method = dist_method)

  diag_dt <- rbindlist(lapply(k_range, function(k) {
    cl <- .cluster_once(z, k = k, method = method, seed = seed + k,
                        nstart = nstart, hclust_method = hclust_method,
                        dist_method = dist_method)
    ## within-cluster sum of squares (to the cluster mean)
    wss <- sum(vapply(sort(unique(cl)), function(c) {
      sub <- z[cl == c, , drop = FALSE]
      if (nrow(sub) < 2) return(0)
      sum(scale(sub, center = TRUE, scale = FALSE)^2)
    }, numeric(1)))
    sil <- NA_real_
    if (requireNamespace("cluster", quietly = TRUE) && length(unique(cl)) > 1) {
      sil <- mean(cluster::silhouette(cl, d_for_sil)[, "sil_width"])
    }
    data.table(k = k, tot_withinss = wss, mean_silhouette = sil)
  }))

  ## distance-to-line elbow
  x <- diag_dt$k; y <- diag_dt$tot_withinss
  x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  dl <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / denom
  elbow <- x[which.max(dl)]
  diag_dt[, elbow_k := elbow]
  setattr(diag_dt, "suggested_k", elbow)

  if (!is.null(file)) {
    subtitle <- NULL
    p <- ggplot(diag_dt, aes(k, tot_withinss))
    if (annotate_k) {
      subtitle <- sprintf("dashed: distance-to-line elbow = k%s", elbow)
      p <- p + geom_vline(xintercept = elbow, linetype = 2, color = "#b2182b")
      if (!is.null(chosen_k)) {
        p <- p + geom_vline(xintercept = chosen_k, linetype = 1,
                            linewidth = 0.9, color = "#2166ac")
        subtitle <- sprintf("%s   |   solid: k used downstream = k%s",
                            subtitle, chosen_k)
      }
    }
    p <- p +
      geom_line(color = "#2b2b2b") + geom_point(size = 2, color = "#2b2b2b") +
      scale_x_continuous(breaks = k_range) +
      theme_bw(base_size = base_size) +
      theme(panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = base_size * 0.85)) +
      labs(x = "k", y = "total within-cluster sum of squares",
           title = "k-means elbow diagnostic", subtitle = subtitle)
    ggsave(file, p, width = 6.2, height = 4)
  }
  diag_dt[]
}

## =============================================================================
## 1. Cluster the data into temporal patterns / modules
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
longitudinal_cluster <- function(mat, meta,
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
                                 n_resample = 50, subsample_frac = 0.8) {

  method <- match.arg(method, c("kmeans", "hierarchical", "pam"))
  meta   <- as.data.table(copy(meta))
  stopifnot(sample_col %in% names(meta), time_col %in% names(meta))
  if (!is.null(group_col)) stopifnot(group_col %in% names(meta))

  ## ---- build the numeric matrix with feature rownames --------------------
  if (!is.null(id_col)) {
    dt  <- as.data.table(mat)
    ids <- as.character(dt[[id_col]])
    keep_cols <- intersect(as.character(meta[[sample_col]]), names(dt))
    m   <- as.matrix(dt[, ..keep_cols])
    rownames(m) <- ids
  } else {
    m <- as.matrix(mat)
    if (is.null(rownames(m)))
      stop("mat has no rownames; supply feature IDs as rownames or via id_col.")
    keep_cols <- intersect(as.character(meta[[sample_col]]), colnames(m))
    m <- m[, keep_cols, drop = FALSE]
  }
  if (length(keep_cols) == 0)
    stop("No overlap between meta[[sample_col]] and the matrix columns.")

  ## keep meta rows matching present samples, in matrix-column order
  meta <- meta[match(keep_cols, meta[[sample_col]])]
  m    <- m[complete.cases(m), , drop = FALSE]
  if (nrow(m) < 2) stop("Fewer than 2 complete-case features to cluster.")

  ## ---- order samples: group (if any) then time then remaining ------------
  ord_cols <- if (is.null(group_col)) c(time_col) else c(group_col, time_col)
  col_order <- do.call(order, c(as.list(meta[, ..ord_cols]),
                                list(seq_len(nrow(meta)))))
  meta <- meta[col_order]
  m    <- m[, meta[[sample_col]], drop = FALSE]

  ## optional replicate collapsing (one column per group x time) ------------
  if (aggregate_replicates) {
    key <- if (is.null(group_col)) as.character(meta[[time_col]])
           else paste(meta[[group_col]], meta[[time_col]], sep = "|")
    key_levels <- unique(key)
    m <- collapse_columns(m, key, key_levels)
    ## rebuild a one-row-per-column meta to match
    meta <- data.table(tmp = key_levels)
    setnames(meta, "tmp", sample_col)
    if (is.null(group_col)) {
      meta[[time_col]] <- as.numeric(key_levels)
    } else {
      parts <- tstrsplit(key_levels, "\\|")
      meta[[group_col]] <- parts[[1]]
      meta[[time_col]]  <- as.numeric(parts[[2]])
    }
  }

  ## ---- z-score ------------------------------------------------------------
  z <- if (is_zscore) m else zscore_rows(m)

  ## ---- cluster ------------------------------------------------------------
  raw <- .cluster_once(z, k = k, method = method, seed = seed,
                       nstart = nstart, hclust_method = hclust_method,
                       dist_method = dist_method)
  k_use <- length(unique(raw))

  ## time grouping for ordering: the numeric time of each column
  time_group  <- as.character(meta[[time_col]])
  time_levels <- as.character(sort(unique(meta[[time_col]])))
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

  ## ---- cluster centroids per (group x) time (for the line plot) ----------
  if (is.null(group_col)) {
    tp_mean  <- collapse_columns(z, as.character(meta[[time_col]]), time_levels)
    centroid <- t(sapply(levels(cluster), function(cl)
      colMeans(tp_mean[cluster == cl, , drop = FALSE], na.rm = TRUE)))
    colnames(centroid) <- time_levels
    centroid_dt <- data.table(Cluster = rownames(centroid), as.data.table(centroid))
  } else {
    gkey <- paste(meta[[group_col]], meta[[time_col]], sep = "|")
    gkey_levels <- unique(gkey)
    tp_mean <- collapse_columns(z, gkey, gkey_levels)
    centroid <- t(sapply(levels(cluster), function(cl)
      colMeans(tp_mean[cluster == cl, , drop = FALSE], na.rm = TRUE)))
    colnames(centroid) <- gkey_levels
    centroid_dt <- data.table(Cluster = rownames(centroid), as.data.table(centroid))
  }

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
                      dist_method = dist_method, hclust_method = hclust_method)
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
    for (f in features) {
      hit <- which(ids == f)
      if (!length(hit)) hit <- grep(f, ids, fixed = TRUE)
      if (length(hit)) { idx <- c(idx, hit); lab <- c(lab, ids[hit]) }
    }
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
## 2. + 5. The ordered heatmap (with metadata bars and feature labels)
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
#'   Defaults to group_col and time_col.
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
#' @param show_stability if TRUE and available, annotate module stability.
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
                                 label_fontsize = 10,
                                 row_title_fontsize = 10,
                                 column_name_fontsize = 9,
                                 annotation_name_fontsize = 10,
                                 legend_fontsize = 10,
                                 title_fontsize = 13,
                                 show_column_names = TRUE,
                                 show_stability = TRUE,
                                 file = NULL, width = 9, height = 13,
                                 title = NULL) {
  stopifnot(inherits(obj, "longi"))
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
  if (is.null(annotation_cols)) {
    annotation_cols <- c(p$group_col, p$time_col)
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
      annotation_name_gp = gpar(fontsize = annotation_name_fontsize),
      annotation_legend_param = .moduleviz_legend_gp(legend_fontsize))
  }

  ## ---- left annotation: module colour + optional stability ---------------
  cluster_cols <- .moduleviz_cluster_cols(levels(cluster), cluster_palette)
  left_args <- list(Cluster = cluster,
                    col = list(Cluster = cluster_cols),
                    show_annotation_name = FALSE,
                    annotation_legend_param = .moduleviz_legend_gp(legend_fontsize),
                    width = unit(4, "mm"))
  if (show_stability && !is.null(obj$stability)) {
    stab_map <- setNames(obj$stability$MeanJaccard, obj$stability$Cluster)
    left_args$Stability <- stab_map[as.character(cluster)]
    left_args$col$Stability <- circlize::colorRamp2(c(0.5, 0.75, 1),
                                                    c("#f7f7f7", "#bdbdbd", "#252525"))
  }
  left_anno <- do.call(rowAnnotation, left_args)

  ## ---- row titles with module size (and stability if present) -------------
  cn <- table(cluster)
  if (!is.null(obj$stability)) {
    stab_map  <- setNames(obj$stability$MeanJaccard, obj$stability$Cluster)
    row_title <- sprintf("%s\n(n=%d, J=%.2f)", levels(cluster),
                         as.integer(cn), stab_map[levels(cluster)])
  } else {
    row_title <- sprintf("%s\n(n=%d)", levels(cluster), as.integer(cn))
  }

  ## ---- right annotation: feature labels -----------------------------------
  lab <- resolve_labels(obj, features = label_features, top_n = top_n_label,
                        rank_table = rank_table, id_field = id_field,
                        value_field = value_field,
                        bigger_is_better = bigger_is_better)
  right_anno <- NULL
  if (nrow(lab)) {
    right_anno <- rowAnnotation(mark = anno_mark(
      at = lab$idx, labels = lab$label,
      labels_gp = gpar(fontsize = label_fontsize, fontface = "italic"),
      link_width = unit(7, "mm")))
  }

  ht <- Heatmap(
    z, name = "z-score", col = col_fun,
    row_split = cluster, cluster_rows = FALSE, cluster_row_slices = FALSE,
    row_title = row_title, row_title_rot = 0,
    row_title_gp = gpar(fontsize = row_title_fontsize),
    row_gap = unit(1.5, "mm"), show_row_names = FALSE,
    column_split = col_split, cluster_columns = FALSE,
    column_gap = unit(1.5, "mm"),
    show_column_names = show_column_names,
    column_names_gp = gpar(fontsize = column_name_fontsize),
    column_title = title,
    column_title_gp = gpar(fontsize = title_fontsize),
    top_annotation = top_anno, left_annotation = left_anno,
    right_annotation = right_anno,
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = c(list(direction = "horizontal",
                                  title_position = "topcenter",
                                  legend_width = unit(4, "cm")),
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
## 3. The per-module pattern / line plot
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
#' @param title plot title; `NULL` drops it entirely.
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
                             base_size = 12,
                             title = "Module patterns",
                             title_fontsize = NULL,
                             xlab = NULL,
                             time_labels = NULL,
                             x_angle = 45,
                             line_width = 1.1,
                             point_size = 2.2) {
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

  if (!is.null(p$group_col)) {
    ## keys look like "Group|Time"
    parts <- tstrsplit(as.character(long$key), "\\|")
    long[, Group := parts[[1]]]
    long[, Time  := as.numeric(parts[[2]])]
    setorder(long, Group, Time)
    long[, Group := factor(Group, levels = unique(Group))]
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
    gg <- ggplot(long, aes(factor(Time), z, group = Group, color = Group)) +
      scale_color_manual(values = group_cols) +
      facet_grid(rows = vars(ClusterF), scales = "free_x")
    n_facet_col <- 1
  } else {
    long[, Time := suppressWarnings(as.numeric(as.character(key)))]
    if (any(is.na(long$Time))) long[, Time := key]      # non-numeric keys
    setorder(long, Cluster, Time)
    cluster_cols <- setNames(
      .moduleviz_discrete_palette(nlevels(obj$cluster), cluster_palette),
      levels(obj$cluster))
    gg <- ggplot(long, aes(factor(Time), z, group = Cluster, color = Cluster)) +
      scale_color_manual(values = cluster_cols, guide = "none") +
      facet_wrap(~ ClusterF, ncol = ncol, strip.position = "right")
    n_facet_col <- ncol
  }

  gg <- gg +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_line(linewidth = line_width) + geom_point(size = point_size) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.text.y.right = element_text(angle = 0),
          plot.title = element_text(size = title_fontsize),
          axis.text.x = element_text(angle = x_angle,
                                     hjust = if (x_angle == 0) 0.5 else 1)) +
    labs(x = if (is.null(xlab)) p$time_col else xlab,
         y = "mean z-score", color = p$group_col, title = title)

  ## optional relabelling of the ordering axis (e.g. 1/2/3 -> naive/mem/exh)
  if (!is.null(time_labels))
    gg <- gg + scale_x_discrete(labels = time_labels)

  if (!is.null(file)) {
    if (is.null(width))  width  <- if (!is.null(p$group_col)) 5.2 else 4.2
    if (is.null(height)) height <- 0.95 * nlevels(obj$cluster) / max(1, n_facet_col) + 0.8
    ggsave(file, gg, width = width, height = height, limitsize = FALSE)
    message("wrote ", file)
  }
  gg
}

## =============================================================================
## 4. Write out memberships / centroids / stability
## =============================================================================

#' Write the cluster memberships (and, optionally, centroids + stability).
#'
#' @param obj a "longi" object.
#' @param file path for the membership table (tab-separated).
#' @param prefix if given (and file is NULL), files are written as
#'   <prefix>_memberships.txt / _centroids.txt / _stability.txt.
#' @param write_extras also write centroids and stability tables.
#' @return invisibly, the membership data.table.
write_memberships <- function(obj, file = NULL, prefix = NULL,
                              write_extras = TRUE) {
  stopifnot(inherits(obj, "longi"))
  if (is.null(file)) {
    if (is.null(prefix)) stop("Supply either `file` or `prefix`.")
    file <- paste0(prefix, "_memberships.txt")
  }
  fwrite(obj$membership, file, sep = "\t", quote = FALSE)
  message("wrote ", file)

  if (write_extras) {
    base <- if (!is.null(prefix)) prefix else sub("\\.[^.]*$", "", file)
    fwrite(obj$centroids, paste0(base, "_centroids.txt"), sep = "\t", quote = FALSE)
    if (!is.null(obj$stability))
      fwrite(obj$stability, paste0(base, "_stability.txt"), sep = "\t", quote = FALSE)
  }
  invisible(obj$membership)
}

## =============================================================================
## 6. Two data sets side by side (RNA + ATAC, RNA + Cut&Run, ...)
## =============================================================================

#' Draw two omics layers side by side, linked by gene.
#'
#' Puts a primary layer (e.g. RNA, one row per gene) next to a secondary layer
#' (e.g. ATAC / Cut&Run peaks) in a single figure, sharing a common set of genes
#' and row order so gene expression and chromatin signal line up per gene.
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
#' @param label_features genes to label down the middle (character vector).
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
dual_omics_heatmap <- function(primary, secondary,
                               primary_name = "RNA", secondary_name = "ATAC",
                               gene_of_primary   = function(id) id,
                               gene_of_secondary = function(id) sub("^.*_", "", id),
                               aggregate = c("mean", "top"),
                               label_features = NULL,
                               col_fun = NULL,
                               zscore_palette = c("rdbu", "viridis"),
                               zscore_limits = c(-2, 2),
                               cluster_palette = "paired",
                               label_fontsize = 10,
                               row_title_fontsize = 10,
                               column_name_fontsize = 9,
                               legend_fontsize = 10,
                               title_fontsize = 13,
                               show_column_names = TRUE,
                               file = NULL, width = 12, height = 13, ...) {
  aggregate <- match.arg(aggregate)
  if (is.null(col_fun))
    col_fun <- .moduleviz_zscore_col_fun(zscore_palette, zscore_limits)

  ## Allow raw matrices: cluster the primary here.
  if (!inherits(primary, "longi")) {
    primary <- longitudinal_cluster(primary, ...)
  }
  p <- primary$params

  ## Secondary z-matrix (either from a longi object, or z-scored raw matrix).
  if (inherits(secondary, "longi")) {
    zsec <- secondary$z
  } else {
    sm <- as.matrix(secondary)
    ## align secondary columns to the primary sample order where possible
    common <- intersect(colnames(primary$z), colnames(sm))
    if (length(common) >= 2) sm <- sm[, common, drop = FALSE]
    zsec <- zscore_rows(sm)
  }

  ## Gene keys.
  gp <- gene_of_primary(rownames(primary$z))
  gs <- gene_of_secondary(rownames(zsec))

  ## Collapse the secondary layer to one row per gene.
  sec_by_gene <- function(g) {
    rows <- which(gs == g)
    if (!length(rows)) return(NULL)
    if (aggregate == "mean") {
      colMeans(zsec[rows, , drop = FALSE], na.rm = TRUE)
    } else {
      strongest <- rows[which.max(apply(abs(zsec[rows, , drop = FALSE]), 1, max))]
      zsec[strongest, ]
    }
  }

  ## Genes present in BOTH layers, kept in the primary's (clustered) row order.
  genes_primary <- gp
  shared        <- genes_primary[genes_primary %in% gs]
  if (!length(shared))
    stop("No genes shared between the two layers; check the gene_of_* mappers.")

  ## Build aligned matrices in primary row order, one row per primary feature
  ## that has a secondary match.
  keep_rows <- which(gp %in% gs)
  z_primary <- primary$z[keep_rows, , drop = FALSE]
  cluster   <- droplevels(primary$cluster[keep_rows])
  g_keep    <- gp[keep_rows]

  z_secondary <- t(vapply(g_keep, function(g) {
    v <- sec_by_gene(g)
    if (is.null(v)) rep(NA_real_, ncol(zsec)) else v
  }, numeric(ncol(zsec))))
  rownames(z_secondary) <- rownames(z_primary)

  ## Column splits (reuse each layer's own sample time ordering).
  make_split <- function(obj_or_mat, zmat) {
    if (inherits(obj_or_mat, "longi") && !is.null(obj_or_mat$params$group_col)) {
      m <- obj_or_mat$meta
      factor(m[[obj_or_mat$params$group_col]], levels = unique(m[[obj_or_mat$params$group_col]]))
    } else NULL
  }
  split_primary   <- make_split(primary, z_primary)
  split_secondary <- if (inherits(secondary, "longi")) make_split(secondary, z_secondary) else NULL

  ## Optional labels down the primary side.
  right_anno <- NULL
  if (!is.null(label_features) && length(label_features)) {
    at <- which(g_keep %in% label_features | rownames(z_primary) %in% label_features)
    if (length(at))
      right_anno <- rowAnnotation(mark = anno_mark(
        at = at, labels = g_keep[at],
        labels_gp = gpar(fontsize = label_fontsize, fontface = "italic"),
        link_width = unit(6, "mm")))
  }

  cluster_cols <- .moduleviz_cluster_cols(levels(cluster), cluster_palette)
  cn <- table(cluster)
  row_title <- sprintf("%s\n(n=%d)", levels(cluster), as.integer(cn))

  ht_primary <- Heatmap(
    z_primary, name = paste0(primary_name, "\nz-score"), col = col_fun,
    row_split = cluster, cluster_rows = FALSE, cluster_row_slices = FALSE,
    row_title = row_title, row_title_rot = 0,
    row_title_gp = gpar(fontsize = row_title_fontsize),
    row_gap = unit(1.5, "mm"), show_row_names = FALSE,
    column_split = split_primary, cluster_columns = FALSE,
    column_title = primary_name, column_gap = unit(1.5, "mm"),
    column_title_gp = gpar(fontsize = title_fontsize),
    show_column_names = show_column_names,
    column_names_gp = gpar(fontsize = column_name_fontsize),
    left_annotation = rowAnnotation(
      Cluster = cluster, col = list(Cluster = cluster_cols),
      show_annotation_name = FALSE,
      annotation_legend_param = .moduleviz_legend_gp(legend_fontsize),
      width = unit(4, "mm")),
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = c(list(direction = "horizontal",
                                  title_position = "topcenter"),
                             .moduleviz_legend_gp(legend_fontsize)))

  ht_secondary <- Heatmap(
    z_secondary, name = paste0(secondary_name, "\nz-score"), col = col_fun,
    cluster_rows = FALSE, show_row_names = FALSE,
    column_split = split_secondary, cluster_columns = FALSE,
    column_title = secondary_name, column_gap = unit(1.5, "mm"),
    column_title_gp = gpar(fontsize = title_fontsize),
    show_column_names = show_column_names,
    column_names_gp = gpar(fontsize = column_name_fontsize),
    right_annotation = right_anno,
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = c(list(direction = "horizontal",
                                  title_position = "topcenter"),
                             .moduleviz_legend_gp(legend_fontsize)))

  combined <- ht_primary + ht_secondary
  draw_it <- function() draw(combined, heatmap_legend_side = "bottom",
                             annotation_legend_side = "bottom", merge_legend = TRUE,
                             column_title = sprintf("%s + %s (linked by gene)",
                                                    primary_name, secondary_name))
  if (!is.null(file)) {
    .moduleviz_write_plot(file, width, height, draw_it)
  } else {
    draw_it()
  }

  invisible(list(primary = primary,
                 z_secondary = z_secondary,
                 genes = g_keep,
                 heatmap = combined))
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
run_longitudinal <- function(mat, meta, out_prefix,
                             label_features = NULL, top_n_label = 5,
                             annotation_cols = NULL, rank_table = NULL,
                             zscore_palette = "rdbu",
                             cluster_palette = "paired",
                             annotation_palette = "set2",
                             heatmap_width = 9, heatmap_height = 13, ...) {
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
