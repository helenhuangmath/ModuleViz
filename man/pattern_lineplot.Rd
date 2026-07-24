\name{pattern_lineplot}
\alias{pattern_lineplot}
\title{Plot mean module trajectories}
\description{
Summarize each module centroid as a faceted line plot across ordered time
points, optionally comparing groups.
}
\usage{
pattern_lineplot(obj, file = NULL, width = NULL, height = NULL, ncol = 1)
}
\arguments{
\item{obj}{A \code{longi} object.}
\item{file}{Optional PDF output path.}
\item{width, height}{Optional plot dimensions.}
\item{ncol}{Facet columns when no group column is present.}
}
\value{A \code{ggplot} object.}
