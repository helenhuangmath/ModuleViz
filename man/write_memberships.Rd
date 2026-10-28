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
