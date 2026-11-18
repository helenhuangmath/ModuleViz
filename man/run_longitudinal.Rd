\name{run_longitudinal}
\alias{run_longitudinal}
\title{Run the full single-omics longitudinal workflow}
\description{
Convenience wrapper that clusters features, draws the heatmap, draws module
trajectory plots, and writes output tables.
}
\usage{
run_longitudinal(mat, meta = NULL, out_prefix, label_features = NULL,
  top_n_label = 5, annotation_cols = NULL, rank_table = NULL,
  zscore_palette = "rdbu", cluster_palette = "paired",
  annotation_palette = "set2", heatmap_width = 5, heatmap_height = 6.5, ...)
}
\arguments{
\item{mat, meta}{Feature matrix and sample metadata.}
\item{out_prefix}{Output file prefix.}
\item{label_features, top_n_label, annotation_cols, rank_table}{Heatmap controls.}
\item{zscore_palette, cluster_palette, annotation_palette}{Palettes forwarded to
\code{longitudinal_heatmap()} and \code{pattern_lineplot()}.}
\item{heatmap_width, heatmap_height}{Heatmap PDF dimensions.}
\item{...}{Passed to \code{longitudinal_cluster()}.}
}
\value{The \code{longi} object invisibly.}
\examples{
data("heatmap_example", package = "ModuleViz")
prefix <- file.path(tempdir(), "moduleviz_example")
obj <- run_longitudinal(heatmap_example$mat, heatmap_example$meta,
                        out_prefix = prefix, time_col = "TimeHr",
                        group_col = "Condition", k = 4, stability = FALSE)
list.files(tempdir(), pattern = "moduleviz_example")
}
