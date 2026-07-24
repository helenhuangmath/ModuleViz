\name{longitudinal_heatmap}
\alias{longitudinal_heatmap}
\alias{resolve_labels}
\title{Draw an ordered longitudinal heatmap}
\description{
Draw a ComplexHeatmap from a \code{longi} object with temporal module splits,
metadata annotations, optional stability bars, and selected feature labels.
}
\usage{
longitudinal_heatmap(obj, annotation_cols = NULL, annotation_colors = NULL,
  label_features = NULL, top_n_label = 0, rank_table = NULL,
  id_field = "ID", value_field = "value", bigger_is_better = TRUE,
  col_fun = NULL, show_stability = TRUE, file = NULL, width = 9,
  height = 13, title = NULL)

resolve_labels(obj, features = NULL, top_n = 0, rank_table = NULL,
  id_field = "ID", value_field = "value", bigger_is_better = TRUE)
}
\arguments{
\item{obj}{A \code{longi} object from \code{longitudinal_cluster()}.}
\item{annotation_cols}{Metadata columns to show as top annotations.}
\item{annotation_colors}{Optional named list of annotation color mappings.}
\item{label_features}{Explicit feature IDs or symbols to label.}
\item{top_n_label, top_n}{Number of automatic labels per cluster.}
\item{rank_table}{Optional feature ranking table.}
\item{id_field, value_field, bigger_is_better}{Ranking-table controls.}
\item{col_fun}{Optional z-score color function.}
\item{show_stability}{Show module stability annotation when available.}
\item{file}{Optional PDF output path.}
\item{width, height}{PDF dimensions in inches.}
\item{title}{Optional heatmap title.}
\item{features}{Explicit features for \code{resolve_labels()}.}
}
\value{
The drawn heatmap object invisibly, or a label-position table for
\code{resolve_labels()}.
}
