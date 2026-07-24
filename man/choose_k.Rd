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
  file = NULL)
}
\arguments{
\item{z}{Row-z-scored matrix.}
\item{k_range}{Integer vector of module counts to test.}
\item{method, nstart, hclust_method, dist_method, seed}{Clustering controls.}
\item{file}{Optional PDF diagnostic plot path.}
}
\value{A data table with diagnostics and attribute \code{suggested_k}.}
