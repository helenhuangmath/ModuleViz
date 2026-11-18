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
  annotation_palette = "set2", label_fontsize = 9,
  row_title_fontsize = 10, column_name_fontsize = 9,
  annotation_name_fontsize = 9, legend_fontsize = 9,
  title_fontsize = 11, show_column_names = TRUE,
  show_stability = FALSE, file = NULL,
  width = 5, height = 6.5, title = NULL, enrichment = NULL,
  enrichment_top_n = 3, enrichment_padj = 0.05, enrichment_fontsize = 9,
  column_labels = c("time", "sample"))

resolve_labels(obj, features = NULL, top_n = 0, rank_table = NULL,
  id_field = "ID", value_field = "value", bigger_is_better = TRUE)
}
\arguments{
\item{enrichment}{Optional module enrichment table from
\code{\link{module_enrichment}} (or any data frame with \code{Cluster},
\code{Term}, and \code{Padj} columns). The top terms of every module are
printed to the right of that module.}
\item{enrichment_top_n}{Number of enriched terms shown per module.}
\item{enrichment_padj}{Adjusted p-value cutoff for a term to be shown.}
\item{enrichment_fontsize}{Font size of the enrichment text.}
\item{obj}{A \code{longi} object from \code{longitudinal_cluster()}.}
\item{annotation_cols}{Metadata columns to show as top annotations. Defaults to
the group column, plus the time column when the time points are not already
printed under the columns.}
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
\item{column_labels}{What to print under each column: \code{"time"} (default)
the time point only, since the group is already shown by the column split and
the top bar; \code{"sample"} the full sample names.}
\item{show_stability}{Show the module stability bar when available (off by
default; the values are in the stability table).}
\item{file}{Optional PDF or PNG output path.}
\item{width, height}{PDF dimensions in inches.}
\item{title}{Optional heatmap title.}
\item{features}{Explicit features for \code{resolve_labels()}.}
}
\value{
The drawn heatmap object invisibly, or a label-position table for
\code{resolve_labels()}.
}
\examples{
data("heatmap_example", package = "ModuleViz")
obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
                            time_col = "TimeHr", group_col = "Condition",
                            k = 4, n_resample = 5)
longitudinal_heatmap(obj, top_n_label = 2)

## print the top enriched gene sets beside each module
genes <- rownames(heatmap_example$mat)
sets <- list(Early_response = genes[1:15], Mid_response = genes[16:30],
             Late_response = genes[31:45])
enr <- module_enrichment(obj, sets, min_size = 3)
longitudinal_heatmap(obj, enrichment = enr)

resolve_labels(obj, features = "Gene001", top_n = 1)
}
