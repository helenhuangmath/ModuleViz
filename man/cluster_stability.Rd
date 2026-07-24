\name{cluster_stability}
\alias{cluster_stability}
\title{Assess cluster stability by row resampling}
\description{
Repeatedly subsample rows, recluster, and estimate each module's mean maximum
Jaccard overlap with resampled clusters.
}
\usage{
cluster_stability(z, cluster, k, method = "kmeans", nstart = 50,
  hclust_method = "ward.D2", dist_method = "euclidean", n_resample = 50,
  subsample_frac = 0.8, seed = 1)
}
\arguments{
\item{z}{Row-z-scored matrix.}
\item{cluster}{Reference cluster labels.}
\item{k, method, nstart, hclust_method, dist_method}{Clustering controls.}
\item{n_resample, subsample_frac, seed}{Resampling controls.}
}
\value{A data table with cluster size and mean Jaccard stability.}
