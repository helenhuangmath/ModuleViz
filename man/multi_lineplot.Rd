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
  ncol = NULL, palette = "okabe_ito", layer_colors = NULL,
  error_bars = TRUE, base_size = 12,
  title = "Module patterns across layers", title_fontsize = NULL,
  xlab = NULL, time_labels = NULL, x_angle = 45, line_width = 1.1,
  point_size = 2.2)
}
\arguments{
\item{x}{A \code{multi_longi} object, for example the value of
\code{multi_heatmap()}.}
\item{file}{Optional output path (any \code{ggsave()} format).}
\item{width, height}{Optional plot dimensions; scaled to the facet grid when
\code{NULL}.}
\item{ncol}{Facet columns when the layers have no group column.}
\item{palette}{Palette for the layer lines. Defaults to \code{"okabe_ito"}.}
\item{layer_colors}{Optional named vector of explicit layer colours, overriding
\code{palette} for the layers it names.}
\item{error_bars}{Draw +/- SD across replicates where there are replicates.}
\item{base_size, title, title_fontsize, xlab, time_labels, x_angle, line_width,
point_size}{As in \code{pattern_lineplot()}.}
}
\details{Layers without a time column in their metadata are left out with a
warning. When layers have a group column, groups become facet columns.}
\value{A ggplot object, also saved when \code{file} is supplied.}
