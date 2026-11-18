\name{multi_lineplot}
\alias{multi_lineplot}
\title{Plot module trajectories across omics layers}
\description{
One facet per module with one line per layer. Each sample contributes its mean
z-score over the module's genes; points are the mean at each time point and
error bars the standard deviation across replicates.
}
\usage{
multi_lineplot(x, file = NULL, width = NULL, height = NULL,
  ncol = 1, split_groups = FALSE, palette = "okabe_ito", layer_colors = NULL,
  error_bars = TRUE, base_size = 10,
  title = NULL, title_fontsize = NULL,
  xlab = NULL, time_labels = NULL, x_angle = 0, line_width = 0.7,
  point_size = 1.5)
}
\arguments{
\item{x}{A \code{multi_longi} object, for example the value of
\code{multi_heatmap()}.}
\item{file}{Optional output path (any \code{ggsave()} format).}
\item{width, height}{Optional plot dimensions; scaled to the facet grid when
\code{NULL}.}
\item{ncol}{Number of facet columns. The default, 1, stacks the modules in a
single column; larger values wrap them into a grid.}
\item{split_groups}{When the layers have a group column: \code{FALSE} (default)
keeps one panel per module and separates the groups by line type;
\code{TRUE} gives every group its own column of panels.}
\item{palette}{Palette for the layer lines. Defaults to \code{"okabe_ito"}.}
\item{layer_colors}{Optional named vector of explicit layer colours, overriding
\code{palette} for the layers it names.}
\item{error_bars}{Draw +/- SD across replicates where there are replicates.}
\item{base_size, title, title_fontsize, xlab, time_labels, x_angle, line_width,
point_size}{As in \code{pattern_lineplot()}.}
}
\details{Layers without a time column in their metadata are left out with a
warning. When layers have a group column, groups are
distinguished by line type, or by facet column with \code{split_groups = TRUE}.}
\value{A ggplot object, also saved when \code{file} is supplied.}
\examples{
data("heatmap_example", package = "ModuleViz")
mat <- heatmap_example$mat
atac <- mat[rep(rownames(mat), each = 2), ]
rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
x <- multi_cluster(list(RNA = mat, ATAC = atac),
                   meta = heatmap_example$meta, sample_col = "SampleID",
                   time_col = "TimeHr", group_col = "Condition",
                   k = 4, stability = FALSE)
multi_lineplot(x)
multi_lineplot(x, split_groups = TRUE)
}
