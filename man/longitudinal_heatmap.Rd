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
  col_fun = NULL, zscore_palette = c("rdbu", "viridis"),
  zscore_limits = c(-2, 2), cluster_palette = "paired",
  annotation_palette = "set2", label_fontsize = 10,
  row_title_fontsize = 10, column_name_fontsize = 9,
  annotation_name_fontsize = 10, legend_fontsize = 10,
  title_fontsize = 13, show_column_names = TRUE,
  show_stability = TRUE, file = NULL,
  width = 9, height = 13, title = NULL)

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
\item{col_fun}{Optional z-score color function. Overrides \code{zscore_palette}.}
\item{zscore_palette}{Continuous palette for the heatmap body: \code{"rdbu"}
(default, blue-white-red) or \code{"viridis"}.}
\item{zscore_limits}{Length-2 z-score range at which the colour scale saturates.}
\item{cluster_palette}{Discrete module color palette. Defaults to \code{"paired"};
\code{"spectral"}, \code{"set2"}, \code{"set3"}, \code{"okabe_ito"},
\code{"nature"}, and \code{"science"} are also available.}
\item{annotation_palette}{Discrete palette for categorical metadata bars such as
condition or genotype. Defaults to \code{"set2"}; use \code{"set3"} for more
than eight levels.}
\item{label_fontsize, row_title_fontsize, column_name_fontsize, annotation_name_fontsize}{Text sizes for heatmap labels and annotations.}
\item{legend_fontsize, title_fontsize}{Text sizes for all legends (titles and
labels) and for the heatmap title.}
\item{show_column_names}{Show sample names along the bottom. Set \code{FALSE}
when sample IDs are long or uninformative and the annotation bars already
identify each column.}
\item{show_stability}{Show module stability annotation when available.}
\item{file}{Optional PDF or PNG output path.}
\item{width, height}{PDF dimensions in inches.}
\item{title}{Optional heatmap title.}
\item{features}{Explicit features for \code{resolve_labels()}.}
}
\value{
The drawn heatmap object invisibly, or a label-position table for
\code{resolve_labels()}.
}
