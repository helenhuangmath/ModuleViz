\name{run_longitudinal}
\alias{run_longitudinal}
\title{Run the full single-omics longitudinal workflow}
\description{
Convenience wrapper that clusters features, draws the heatmap, draws module
trajectory plots, and writes output tables.
}
\usage{
run_longitudinal(mat, meta, out_prefix, label_features = NULL,
  top_n_label = 5, annotation_cols = NULL, rank_table = NULL,
  heatmap_width = 9, heatmap_height = 13, ...)
}
\arguments{
\item{mat, meta}{Feature matrix and sample metadata.}
\item{out_prefix}{Output file prefix.}
\item{label_features, top_n_label, annotation_cols, rank_table}{Heatmap controls.}
\item{heatmap_width, heatmap_height}{Heatmap PDF dimensions.}
\item{...}{Passed to \code{longitudinal_cluster()}.}
}
\value{The \code{longi} object invisibly.}
