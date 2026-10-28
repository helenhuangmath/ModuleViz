\name{multi_cluster}
\alias{multi_cluster}
\alias{print.multi_longi}
\title{Cluster several omics layers into shared modules, linked by gene}
\description{
Align any number of omics layers (for example RNA, ATAC, H3K27ac, and
H3K27me3) by gene and assign every gene to one temporal module. Each layer is
z-scored within itself and collapsed to one row per gene, so layers with
different feature types (genes, peaks, CUT&RUN regions) line up row for row.
}
\usage{
multi_cluster(layers, meta = NULL, sample_col = "SampleID",
  time_col = "TimeHr", group_col = NULL, cluster_by = 1, gene_of = NULL,
  aggregate = c("mean", "top"), gene_set = c("intersect", "reference"),
  balance_layers = TRUE, ...)
}
\arguments{
\item{layers}{Named list of at least two layers. Each entry is a \code{longi}
object or a raw features x samples matrix. Names become panel titles.}
\item{meta}{Sample metadata for raw-matrix layers: one data frame shared by all
layers, or a named list with one data frame per layer. \code{longi} layers carry
their own metadata.}
\item{sample_col, time_col, group_col}{Metadata columns for raw-matrix layers,
as in \code{longitudinal_cluster()}.}
\item{cluster_by}{Which layer defines the modules: a layer name or index
(default: the first layer), or \code{"joint"} to cluster the concatenated
per-gene profiles of every layer together.}
\item{gene_of}{Function mapping row IDs to gene symbols, either one function for
every layer or a list with one per layer (by name or position). By default the
first layer's IDs are used as-is and every other layer takes the text after the
last underscore, so \code{peak1243_Pdcd1} maps to \code{Pdcd1}.}
\item{aggregate}{How to collapse several features per gene: \code{"mean"}
averages them, \code{"top"} keeps the feature with the strongest |z|.}
\item{gene_set}{\code{"intersect"} keeps genes present in every layer;
\code{"reference"} keeps every row of the \code{cluster_by} layer and leaves
layers without that gene as \code{NA}.}
\item{balance_layers}{Joint clustering only: weight each layer's columns by
1/sqrt(number of columns) so a layer with many samples does not dominate the
distance.}
\item{...}{Passed to \code{longitudinal_cluster()} (\code{k}, \code{method},
\code{stability}, ...) when the clustering layer is a raw matrix, and always
for \code{cluster_by = "joint"}.}
}
\details{
With \code{cluster_by} set to a layer, that layer alone defines the modules and
row order and the other layers are displayed in the same order, as in
\code{dual_heatmap()}. With \code{cluster_by = "joint"} a module groups
genes whose trajectories agree across all layers. Joint clustering needs a time
column in every layer's metadata. To collapse replicates, pass per-layer
\code{longitudinal_cluster(..., aggregate_replicates = TRUE)} objects.
}
\value{An object of class \code{multi_longi}: a list with \code{layers}
(per-layer aligned z-score matrix and metadata), \code{cluster}, \code{genes},
\code{membership}, \code{stability}, \code{reference} (the \code{longi} object
the modules came from), and \code{params}.}
\examples{
data("heatmap_example", package = "ModuleViz")
mat <- heatmap_example$mat
atac <- mat[rep(rownames(mat), each = 2), ]
rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
x <- multi_cluster(list(RNA = mat, ATAC = atac),
                         meta = heatmap_example$meta, sample_col = "SampleID",
                         time_col = "TimeHr", group_col = "Condition",
                         k = 4, stability = FALSE)
x
}
