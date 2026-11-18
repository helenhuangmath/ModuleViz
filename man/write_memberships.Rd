\name{write_memberships}
\alias{write_memberships}
\title{Write module membership and summary tables}
\description{
Export module membership, centroid, and stability tables from a \code{longi} or \code{multi_longi}
object.
}
\usage{
write_memberships(obj, file = NULL, prefix = NULL, write_extras = TRUE)
}
\arguments{
\item{obj}{A \code{longi} or \code{multi_longi} object. A \code{multi_longi} has no centroid table, so only membership and stability are written.}
\item{file}{Optional membership table path.}
\item{prefix}{Output prefix used when \code{file} is not supplied.}
\item{write_extras}{Also write centroids and stability tables.}
}
\value{The membership table invisibly.}
\examples{
data("heatmap_example", package = "ModuleViz")
obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
                            time_col = "TimeHr", k = 4, n_resample = 5)
prefix <- file.path(tempdir(), "moduleviz_tables")
write_memberships(obj, prefix = prefix)
list.files(tempdir(), pattern = "moduleviz_tables")
}
