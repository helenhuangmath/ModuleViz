\name{dual_heatmap}
\alias{dual_heatmap}
\title{Draw paired omics heatmaps aligned by gene}
\description{
The two-layer case of \code{multi_heatmap()}: place a primary clustered
omics layer next to a secondary layer collapsed and
aligned by gene, so paired trajectories can be compared in one figure.
}
\usage{
dual_heatmap(primary, secondary, primary_name = "RNA",
  secondary_name = "ATAC", gene_of_primary = function(id) id,
  gene_of_secondary = function(id) sub("^.*_", "", id),
  aggregate = c("mean", "top"), label_features = NULL, col_fun = NULL,
  zscore_palette = c("rdbu", "viridis"), zscore_limits = c(-2, 2),
  cluster_palette = "paired", label_fontsize = 15,
  row_title_fontsize = 16, column_name_fontsize = 13,
  legend_fontsize = 14, title_fontsize = 17, show_column_names = TRUE,
  file = NULL, width = 10.5, height = 13, ...)
}
\arguments{
\item{primary, secondary}{\code{longi} objects or raw matrices.}
\item{primary_name, secondary_name}{Layer names used in titles and legends.}
\item{gene_of_primary, gene_of_secondary}{Functions mapping row IDs to gene symbols.}
\item{aggregate}{How to collapse multiple secondary features per gene.}
\item{label_features}{Optional genes to label.}
\item{col_fun}{Optional z-score color function. Overrides \code{zscore_palette}.}
\item{zscore_palette}{Continuous palette shared by both layers: \code{"rdbu"}
(default) or \code{"viridis"}.}
\item{zscore_limits}{Length-2 z-score range at which the colour scale saturates.}
\item{cluster_palette}{Discrete module color palette. Defaults to \code{"paired"}.}
\item{label_fontsize, row_title_fontsize, column_name_fontsize}{Text sizes for labels and heatmap titles.}
\item{legend_fontsize, title_fontsize}{Text sizes for the legends and the
per-layer column titles.}
\item{show_column_names}{Show sample names along the bottom of both layers.}
\item{file}{Optional PDF or PNG output path.}
\item{width, height}{PDF dimensions in inches.}
\item{...}{Passed to \code{longitudinal_cluster()} when primary is raw.}
}
\value{A list containing the primary object, aligned secondary matrix, genes, and heatmap.}
