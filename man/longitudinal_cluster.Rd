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
longitudinal_cluster(mat, meta = NULL, sample_col = "SampleID",
  time_col = "TimeHr",
  group_col = NULL, id_col = NULL, k = 8, method = "kmeans",
  is_zscore = FALSE, aggregate_replicates = FALSE, seed = 1, nstart = 50,
  hclust_method = "ward.D2", dist_method = "euclidean", stability = TRUE,
  n_resample = 50, subsample_frac = 0.8, assay = NULL)

zscore_rows(mat)

collapse_columns(mat, grp, group_levels = NULL)
}
\arguments{
\item{assay}{Assay to use when \code{mat} is a \code{SummarizedExperiment}.
Defaults to the first of \code{logcpm}, \code{logcounts}, \code{normcounts},
\code{vst}, \code{rlog} that is present, otherwise the first assay.}
\item{mat}{Features-by-samples data: a numeric matrix or data frame, or a
\code{SummarizedExperiment} (including \code{SingleCellExperiment}, e.g. the
result of \code{\link{pseudobulk_cells}}).}
\item{meta}{Sample metadata with one row per sample. Optional for a
\code{SummarizedExperiment}, whose \code{colData} is used by default.}
\item{sample_col}{Metadata column containing sample IDs matching the column
names of \code{mat}. For a \code{SummarizedExperiment} without \code{meta} the
column names are used directly.}
\item{time_col}{Metadata column giving the time point of each sample. Numeric
time is ordered numerically, a factor by its levels, and character values by
first appearance.}
\item{group_col}{Optional metadata column for treatment, genotype, or condition.}
\item{id_col}{Optional feature-ID column when \code{mat} is a data frame.}
\item{k}{Number of modules.}
\item{method}{Clustering algorithm: \code{"kmeans"} (default),
\code{"hierarchical"}, \code{"pam"}, \code{"clara"}, or \code{"fuzzy"}. See
Details.}
\item{is_zscore}{Set \code{TRUE} when values are already row z-scores.}
\item{aggregate_replicates}{Average replicate columns within group-by-time
before clustering. Metadata columns that are constant within each group-by-time
cell are kept.}
\item{seed, nstart, hclust_method, dist_method}{Clustering controls.
\code{seed} makes the randomised algorithms reproducible without altering the
global random-number stream; \code{NULL} uses the current RNG state.
\code{nstart} applies to k-means, \code{hclust_method} to hierarchical
clustering, and \code{dist_method} (\code{"euclidean"}, \code{"correlation"},
or any \code{\link[stats]{dist}} method) to hierarchical, PAM, and fuzzy
clustering.}
\item{stability, n_resample, subsample_frac}{Cluster stability controls.}
\item{grp}{Column grouping vector for \code{collapse_columns()}.}
\item{group_levels}{Optional output group order.}
}
\details{
The clustering algorithms differ in what they assume and how they scale:
\describe{
\item{\code{"kmeans"}}{Fast and the usual choice for z-scored trajectories.}
\item{\code{"hierarchical"}}{Deterministic agglomerative clustering
(\code{hclust_method}, \code{dist_method}); useful with correlation distance.}
\item{\code{"pam"}}{Partitioning around medoids; more robust to outlying
features than k-means, but needs a full distance matrix.}
\item{\code{"clara"}}{PAM on repeated sub-samples; scales to tens of thousands
of features.}
\item{\code{"fuzzy"}}{Fuzzy c-means (\code{\link[cluster]{fanny}}). Features
are assigned to their strongest module and the membership table gains a
\code{Membership} column (0-1) showing how firmly each feature belongs there.}
}
}
\value{
A \code{longi} object containing ordered z-scores, cluster labels, membership
tables, centroids, stability, metadata, and parameters.
}
\seealso{\code{\link{choose_k}}, \code{\link{longitudinal_heatmap}},
\code{\link{pattern_lineplot}}, \code{\link{pseudobulk_cells}}}
\examples{
data("heatmap_example", package = "ModuleViz")
obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
                            sample_col = "SampleID", time_col = "TimeHr",
                            group_col = "Condition", k = 4, n_resample = 5)
obj
head(obj$membership)

## another algorithm: fuzzy c-means adds a membership score
fz <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
                           time_col = "TimeHr", group_col = "Condition",
                           k = 4, method = "fuzzy", stability = FALSE)
head(fz$membership)

## a SummarizedExperiment works directly
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = heatmap_example$mat),
    colData = heatmap_example$meta)
colnames(se) <- heatmap_example$meta$SampleID
longitudinal_cluster(se, time_col = "TimeHr", k = 4, stability = FALSE)

z <- zscore_rows(heatmap_example$mat)
collapse_columns(z, heatmap_example$meta$TimeHr)[1:3, ]
}
