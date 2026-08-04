\name{choose_k}
\alias{choose_k}
\title{Choose a module count with elbow and silhouette diagnostics}
\description{
Scan candidate module counts and report total within-cluster sum of squares,
mean silhouette width, and a simple distance-to-line elbow suggestion.
}
\usage{
choose_k(z, k_range = 2:20, method = "kmeans", nstart = 25,
  hclust_method = "ward.D2", dist_method = "euclidean", seed = 100,
  file = NULL, chosen_k = NULL, annotate_k = TRUE, base_size = 11)
}
\arguments{
\item{z}{Row-z-scored matrix.}
\item{k_range}{Integer vector of module counts to test.}
\item{method, nstart, hclust_method, dist_method, seed}{Clustering controls.}
\item{file}{Optional PDF diagnostic plot path.}
\item{chosen_k}{Optional module count actually used downstream. When supplied it
is drawn as a solid line beside the dashed elbow, keeping the diagnostic figure
consistent with the module count used in the rest of the analysis. The elbow is
only a heuristic and is frequently not the \code{k} worth reporting.}
\item{annotate_k}{If \code{FALSE}, draw the bare elbow curve with no marker
lines and no subtitle.}
\item{base_size}{Base font size for the diagnostic plot.}
}
\value{A data table with diagnostics and attribute \code{suggested_k}.}
\note{
The elbow is a heuristic, not an answer. Prefer the \code{k} that separates
programmes you can interpret: a \code{k} that merges two biologically distinct
modules is worse than one that splits a module in two, because a merged module
reports the average of programmes that do not agree.
}
