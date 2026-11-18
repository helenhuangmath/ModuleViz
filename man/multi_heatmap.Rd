\name{multi_heatmap}
\alias{multi_heatmap}
\title{Draw several omics layers side by side, linked by gene}
\description{
One heatmap panel per layer (for example RNA | ATAC | H3K27ac | H3K27me3), all
sharing the module split and row order from \code{multi_cluster()}. The
layers share one z-score colour scale, so one legend covers every panel; each
panel has its own metadata annotation bars, coloured consistently across
layers.
}
\usage{
multi_heatmap(x, annotation_cols = NULL, annotation_colors = NULL,
  label_features = NULL, top_n_label = 0, col_fun = NULL,
  zscore_palette = c("rdbu", "viridis"), zscore_limits = c(-2, 2),
  cluster_palette = "paired", annotation_palette = "set2",
  na_col = "grey85", label_fontsize = 9, row_title_fontsize = 10,
  column_name_fontsize = 9, annotation_name_fontsize = 9,
  legend_fontsize = 9, title_fontsize = 11, show_column_names = TRUE,
  layer_gap = 5, title = NULL, file = NULL, width = NULL, height = 6.5, ...,
  enrichment = NULL, enrichment_top_n = 3, enrichment_padj = 0.05,
  enrichment_fontsize = 9, column_labels = c("time", "sample"),
  layer_widths = c("equal", "samples"))
}
\arguments{
\item{enrichment}{Optional module enrichment table from
\code{\link{module_enrichment}} (or any data frame with \code{Cluster},
\code{Term}, and \code{Padj} columns). The top terms of every module are
printed to the right of the last panel.}
\item{enrichment_top_n}{Number of enriched terms shown per module.}
\item{enrichment_padj}{Adjusted p-value cutoff for a term to be shown.}
\item{enrichment_fontsize}{Font size of the enrichment text.}
\item{x}{A \code{multi_longi} object, or a named list of layers, which is passed
to \code{multi_cluster()} together with \code{...}.}
\item{annotation_cols}{Metadata columns shown as top bars on each layer.
Defaults to each layer's group column (plus the time column when the time
points are not printed under the columns); \code{character(0)} for none.}
\item{annotation_colors}{Optional named list mapping annotation to colours.}
\item{label_features}{Genes or IDs to label on the right of the last panel.}
\item{top_n_label}{Automatic labels per module, ranked by the largest |z| across
all layers (0 = none).}
\item{col_fun}{Optional z-score colour function. Overrides \code{zscore_palette}.}
\item{zscore_palette}{\code{"rdbu"} (default) or \code{"viridis"}, shared by every
layer.}
\item{zscore_limits}{Length-2 z-score range at which the colour scale saturates.}
\item{cluster_palette}{Discrete module colour palette.}
\item{annotation_palette}{Discrete palette for categorical metadata bars.}
\item{na_col}{Colour for genes missing from a layer (\code{gene_set = "reference"}).}
\item{label_fontsize, row_title_fontsize, column_name_fontsize}{Text sizes for
gene labels, module names, and sample names.}
\item{annotation_name_fontsize, legend_fontsize, title_fontsize}{Text sizes for
annotation names, legends, and panel / figure titles.}
\item{show_column_names}{Show column labels along the bottom of every layer.}
\item{column_labels}{\code{"time"} (default) prints only the time point under
each column; \code{"sample"} prints the full sample names.}
\item{layer_gap}{Gap between the layer panels, in mm.}
\item{layer_widths}{\code{"equal"} (default) gives every layer the same panel
width; \code{"samples"} makes each panel as wide as its number of samples.}
\item{title}{Optional overall title; none by default, since each panel is
already headed by its layer name.}
\item{file}{Optional PDF or PNG output path.}
\item{width, height}{Size in inches; \code{width} defaults to 2.5 + 2 per layer.}
\item{...}{Passed to \code{multi_cluster()} when \code{x} is a list.}
}
\value{The \code{multi_longi} object invisibly, with the drawn
\code{HeatmapList} in \code{$heatmap}.}
\examples{
data("heatmap_example", package = "ModuleViz")
mat <- heatmap_example$mat
atac <- mat[rep(rownames(mat), each = 2), ]
rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
k27 <- -mat
rownames(k27) <- paste0("region", seq_len(nrow(k27)), "_", rownames(k27))
res <- multi_heatmap(list(RNA = mat, ATAC = atac, H3K27me3 = k27),
                           meta = heatmap_example$meta,
                           sample_col = "SampleID", time_col = "TimeHr",
                           group_col = "Condition", k = 4, stability = FALSE,
                           file = tempfile(fileext = ".pdf"))
}
