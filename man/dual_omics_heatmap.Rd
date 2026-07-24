\name{dual_omics_heatmap}
\alias{dual_omics_heatmap}
\title{Draw paired omics heatmaps aligned by gene}
\description{
Place a primary clustered omics layer next to a secondary layer collapsed and
aligned by gene, so paired trajectories can be compared in one figure.
}
\usage{
dual_omics_heatmap(primary, secondary, primary_name = "RNA",
  secondary_name = "ATAC", gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = c("mean", "top"), label_features = NULL, col_fun = NULL,
  file = NULL, width = 12, height = 13, ...)
}
\arguments{
\item{primary, secondary}{\code{longi} objects or raw matrices.}
\item{primary_name, secondary_name}{Layer names used in titles and legends.}
\item{gene_of_primary, gene_of_secondary}{Functions mapping row IDs to gene symbols.}
\item{aggregate}{How to collapse multiple secondary features per gene.}
\item{label_features}{Optional genes to label.}
\item{col_fun}{Optional z-score color function.}
\item{file}{Optional PDF output path.}
\item{width, height}{PDF dimensions in inches.}
\item{...}{Passed to \code{longitudinal_cluster()} when primary is raw.}
}
\value{A list containing the primary object, aligned secondary matrix, genes, and heatmap.}
