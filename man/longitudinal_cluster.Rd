\name{longitudinal_cluster}
\alias{longitudinal_cluster}
\alias{zscore_rows}
\alias{collapse_columns}
\title{Cluster longitudinal features into temporal modules}
\description{
Z-score feature trajectories, cluster them into modules, order modules by their
temporal peak, and return a \code{longi} object for plotting and export.
}
\usage{
longitudinal_cluster(mat, meta, sample_col = "SampleID", time_col = "TimeHr",
  group_col = NULL, id_col = NULL, k = 8, method = "kmeans",
  is_zscore = FALSE, aggregate_replicates = FALSE, seed = 1, nstart = 50,
  hclust_method = "ward.D2", dist_method = "euclidean", stability = TRUE,
  n_resample = 50, subsample_frac = 0.8)

zscore_rows(mat)

collapse_columns(mat, grp, group_levels = NULL)
}
\arguments{
\item{mat}{Numeric matrix or data frame with features in rows and samples in columns.}
\item{meta}{Sample metadata with one row per sample.}
\item{sample_col}{Metadata column containing sample IDs.}
\item{time_col}{Numeric metadata column containing ordered time values.}
\item{group_col}{Optional metadata column for treatment, genotype, or condition.}
\item{id_col}{Optional feature-ID column when \code{mat} is a data frame.}
\item{k}{Number of modules.}
\item{method}{One of \code{"kmeans"}, \code{"hierarchical"}, or \code{"pam"}.}
\item{is_zscore}{Set \code{TRUE} when values are already row z-scores.}
\item{aggregate_replicates}{Average replicate columns within group-by-time before clustering.}
\item{seed, nstart, hclust_method, dist_method}{Clustering controls.}
\item{stability, n_resample, subsample_frac}{Cluster stability controls.}
\item{grp}{Column grouping vector for \code{collapse_columns()}.}
\item{group_levels}{Optional output group order.}
}
\value{
A \code{longi} object containing ordered z-scores, cluster labels, membership
tables, centroids, stability, metadata, and parameters.
}
