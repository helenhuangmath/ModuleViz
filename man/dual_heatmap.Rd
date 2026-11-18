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
  aggregate = c("mean", "top"), label_features = NULL,
  col_fun = NULL,
  zscore_palette = c("rdbu", "viridis"), zscore_limits = c(-2, 2),
  cluster_palette = "paired", label_fontsize = 9,
  row_title_fontsize = 10, column_name_fontsize = 9,
  legend_fontsize = 9, title_fontsize = 11, show_column_names = TRUE,
  file = NULL, width = 6.5, height = 6.5, ..., enrichment = NULL,
  enrichment_top_n = 3, enrichment_padj = 0.05, annotation_cols = NULL,
  annotation_colors = NULL)
}
\arguments{
\item{enrichment, enrichment_top_n, enrichment_padj}{Optional module enrichment
table printed beside the modules; see \code{\link{multi_heatmap}}.}
\item{primary, secondary}{\code{longi} objects or raw matrices.}
\item{primary_name, secondary_name}{Layer names used in titles and legends.}
\item{gene_of_primary, gene_of_secondary}{Functions mapping row IDs to gene symbols.}
\item{aggregate}{How to collapse multiple secondary features per gene.}
\item{label_features}{Optional genes to label.}
\item{annotation_cols, annotation_colors}{Top annotation bars and their colours,
as in \code{multi_heatmap()}; by default the group column of each layer that
carries sample metadata.}
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
\examples{
data("heatmap_example", package = "ModuleViz")
rna <- heatmap_example$mat
atac <- rna[rep(rownames(rna), each = 2), ]
rownames(atac) <- paste0("peak", seq_len(nrow(atac)), "_", rownames(atac))
obj <- longitudinal_cluster(rna, heatmap_example$meta, time_col = "TimeHr",
                            group_col = "Condition", k = 4, stability = FALSE)
res <- dual_heatmap(obj, atac, label_features = c("Gene001", "Gene031"))
names(res)
}
