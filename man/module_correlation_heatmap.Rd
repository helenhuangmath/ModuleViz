\name{module_correlation_heatmap}
\alias{module_correlation_heatmap}
\alias{module_correlation_summary}
\title{Gene-gene correlation heatmaps per temporal module}
\description{
Correlate the features within each temporal module against each other across
the ordered samples and draw the resulting correlation matrix as a heatmap.
A coherent module appears as a solid block of high correlation; visible
sub-blocks mean the module is averaging over more than one programme.
\code{module_correlation_summary()} returns the mean within-module correlation
for every module without drawing anything.
}
\usage{
module_correlation_heatmap(obj, clusters = NULL,
  method = c("pearson", "spearman"), max_genes = Inf,
  select = c("top", "random"), cluster_genes = TRUE,
  corr_palette = "rdbu", corr_limits = c(-1, 1), col_fun = NULL,
  show_gene_names = NULL, label_top_n = 30, label_fontsize = 8,
  legend_fontsize = 10, title_fontsize = 11, file = NULL,
  width = 7, height = 6.5, title = NULL, seed = 1)

module_correlation_summary(obj, method = c("pearson", "spearman"),
  max_genes = Inf, select = c("top", "random"), seed = 1)
}
\arguments{
\item{obj}{A \code{longi} object from \code{longitudinal_cluster()}.}
\item{clusters}{Module names to draw. Defaults to all modules.}
\item{method}{Correlation method, \code{"pearson"} or \code{"spearman"}.}
\item{max_genes}{Optional cap on features per module. Defaults to \code{Inf},
so every feature in the module is used; set a finite value to subsample large
modules for speed.}
\item{select}{Subsampling rule: \code{"top"} keeps the strongest-signal
features (largest max |z|), \code{"random"} draws an unbiased sample.}
\item{cluster_genes}{Order features by hierarchical clustering on correlation
distance, which makes sub-blocks visible.}
\item{corr_palette}{Continuous palette: \code{"rdbu"} (default) or
\code{"viridis"}.}
\item{corr_limits}{Length-2 correlation range for the colour scale.}
\item{col_fun}{Optional explicit colour function; overrides \code{corr_palette}.}
\item{show_gene_names}{Show every feature name. Defaults to \code{TRUE} when 60
or fewer features are shown, and \code{FALSE} beyond that, where printing all
names would be unreadable.}
\item{label_top_n}{When names are too crowded to print in full, mark this many
features instead, chosen as those with the largest \eqn{|z|} - the features
moving most across the series. Set \code{0} to disable marking.}
\item{label_fontsize}{Font size for feature names.}
\item{legend_fontsize, title_fontsize}{Font sizes for the legend and title.}
\item{file}{Optional output path. A \code{.pdf} path collects all modules into
one multi-page PDF; any other extension writes one file per module with the
module name inserted before the extension.}
\item{width, height}{Figure size in inches.}
\item{title}{Optional title prefix; the module name, size, and mean correlation
are appended.}
\item{seed}{RNG seed used when \code{select = "random"}.}
}
\value{
For \code{module_correlation_heatmap()}, invisibly a named list with one entry
per module holding \code{cor}, \code{mean_cor}, \code{n_total}, \code{n_shown},
\code{labelled} (the feature names actually drawn), and \code{heatmap}. For \code{module_correlation_summary()}, a
\code{data.table} with columns \code{Cluster}, \code{N}, \code{NShown}, and
\code{MeanCor}.
}
\examples{
data("heatmap_example", package = "ModuleViz")
obj <- longitudinal_cluster(
  heatmap_example$mat, heatmap_example$meta,
  sample_col = "SampleID", time_col = "TimeHr",
  k = 4, stability = FALSE
)
module_correlation_summary(obj)
}
