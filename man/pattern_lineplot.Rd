\name{pattern_lineplot}
\alias{pattern_lineplot}
\title{Plot mean module trajectories}
\description{
Summarize each module centroid as a faceted line plot across ordered time
points, optionally comparing groups.
}
\usage{
pattern_lineplot(obj, file = NULL, width = NULL, height = NULL, ncol = 1,
  palette = "set2", cluster_palette = "paired", group_colors = NULL,
  base_size = 12,
  title = "Module patterns", title_fontsize = NULL, xlab = NULL,
  time_labels = NULL, x_angle = 45, line_width = 1.1, point_size = 2.2)
}
\arguments{
\item{obj}{A \code{longi} object.}
\item{file}{Optional PDF output path.}
\item{width, height}{Optional plot dimensions.}
\item{ncol}{Facet columns when no group column is present.}
\item{palette}{Palette for the group / condition lines. Defaults to
\code{"set2"}; \code{"set3"}, \code{"paired"}, \code{"nature"}, and
\code{"science"} are also available.}
\item{cluster_palette}{Palette for the per-module lines when the object has no
group column. Defaults to \code{"paired"}.}
\item{group_colors}{Optional named character vector of explicit colours for the
group lines, e.g. \code{c(Control = "grey70", Stimulated = "firebrick")}.
Overrides \code{palette} for the levels it names; unnamed levels keep their
palette colour.}
\item{base_size}{Base ggplot font size.}
\item{title}{Plot title; \code{NULL} drops it entirely.}
\item{title_fontsize}{Title size in points. Defaults to \code{base_size}, so the
title stays proportionate rather than inflating with the base size.}
\item{xlab}{X-axis label; defaults to the name of \code{time_col}.}
\item{time_labels}{Optional relabelling of the x-axis ticks, for example
\code{c("1" = "Naive", "2" = "Memory", "3" = "Exhausted")} when the ordering
column is an integer rank rather than real time.}
\item{x_angle}{Rotation of the x tick labels in degrees. Use \code{0} when the
labels are short enough to sit flat.}
\item{line_width, point_size}{Line and point size.}
}
\value{A \code{ggplot} object.}
