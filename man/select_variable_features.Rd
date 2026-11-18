\name{select_variable_features}
\alias{select_variable_features}
\title{Keep the most variable features}
\description{
Restrict a matrix or \code{SummarizedExperiment} to the features with the
largest variance across samples. Most features of a transcriptome-wide dataset
are flat across a time course, so this is the usual step between
\code{\link{pseudobulk_cells}} (or a bulk matrix) and
\code{\link{longitudinal_cluster}} when no differential-expression result is
available to choose features with.
}
\usage{
select_variable_features(x, n = 2000, assay = NULL)
}
\arguments{
\item{x}{A features-by-samples matrix or a \code{SummarizedExperiment}.}
\item{n}{Number of features to keep.}
\item{assay}{Assay used to rank the features when \code{x} is a
\code{SummarizedExperiment}; chosen as in \code{longitudinal_cluster()} by
default.}
}
\value{\code{x} restricted to its \code{n} most variable features, most variable
first.}
\examples{
data("heatmap_example", package = "ModuleViz")
top <- select_variable_features(heatmap_example$mat, n = 20)
dim(top)
}
