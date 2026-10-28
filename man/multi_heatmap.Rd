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
  na_col = "grey85", label_fontsize = 15, row_title_fontsize = 16,
  column_name_fontsize = 13, annotation_name_fontsize = 14,
  legend_fontsize = 14, title_fontsize = 17, show_column_names = TRUE,
  layer_gap = 10, title = NULL, file = NULL, width = NULL,
  height = 13, ...)
}
\arguments{
\item{x}{A \code{multi_longi} object, or a named list of layers, which is passed
to \code{multi_cluster()} together with \code{...}.}
\item{annotation_cols}{Metadata columns shown as top bars on each layer.
Defaults to each layer's group and time columns; \code{character(0)} for none.}
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
\item{show_column_names}{Show sample names along the bottom of every layer.}
\item{layer_gap}{Gap between the layer panels, in mm.}
\item{title}{Overall title. \code{NULL} gives "A + B + C (linked by gene)";
\code{FALSE} drops it.}
\item{file}{Optional PDF or PNG output path.}
\item{width, height}{Size in inches; \code{width} defaults to 3.5 + 3.5 per layer.}
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
