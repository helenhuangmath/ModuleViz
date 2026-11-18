\name{module_enrichment}
\alias{module_enrichment}
\alias{enrichment_dotplot}
\title{Gene-set enrichment of temporal modules}
\description{
\code{module_enrichment()} tests every module for over-representation of each
gene set with a one-sided hypergeometric test. The result can be passed to the
\code{enrichment} argument of \code{\link{longitudinal_heatmap}},
\code{\link{multi_heatmap}}, and \code{\link{dual_heatmap}} to print the top
terms beside each module, or drawn on its own with \code{enrichment_dotplot()}.
}
\usage{
module_enrichment(obj, gene_sets, universe = NULL, gene_of = NULL,
  min_size = 5, max_size = 500, min_overlap = 2, p_adjust_method = "BH")

enrichment_dotplot(enrichment, top_n = 5, padj_cutoff = 0.05,
  base_size = 10, title = NULL, file = NULL, width = NULL,
  height = NULL)
}
\arguments{
\item{obj}{A \code{longi} or \code{multi_longi} object.}
\item{gene_sets}{Gene sets to test: a named list of gene vectors, a two-column
data frame (term, gene), or the path to a GMT file.}
\item{universe}{Background genes. Defaults to every gene in \code{obj}, which
is appropriate when the object holds all tested genes. Supply the full set of
measured genes when \code{obj} was built from a pre-selected subset such as
differentially expressed genes.}
\item{gene_of}{Optional function mapping feature IDs to gene symbols, for
objects whose rows are not genes, e.g.
\code{function(id) sub("^.*_", "", id)} for \code{peak12_Cd8a}. Ignored for
\code{multi_longi} objects, which already carry gene symbols.}
\item{min_size, max_size}{Gene sets with fewer or more genes in the universe
are skipped.}
\item{min_overlap}{Module / gene-set pairs sharing fewer genes are not tested.}
\item{p_adjust_method}{Multiple-testing correction passed to
\code{\link[stats]{p.adjust}}, applied across all tested pairs.}
\item{enrichment}{The table returned by \code{module_enrichment()}, or any data
frame with \code{Cluster}, \code{Term}, \code{Padj}, and
\code{FoldEnrichment} columns (for example, results converted from another
enrichment tool).}
\item{top_n}{Number of terms shown per module, ranked by adjusted p-value.}
\item{padj_cutoff}{Terms with a larger adjusted p-value are not shown.}
\item{base_size}{ggplot base font size.}
\item{title}{Optional plot title; none by default.}
\item{file}{Optional output path (any \code{ggsave()} format).}
\item{width, height}{Plot size in inches; scaled to the content when
\code{NULL}.}
}
\details{
The test needs no annotation package: any collection of gene sets works,
including pathway databases exported as GMT, transcription-factor targets, or
custom marker lists, for any organism. Gene identifiers in \code{gene_sets}
must use the same naming as the features of \code{obj}.

The choice of \code{universe} matters. Testing modules of pre-selected genes
against only those genes asks "is this set concentrated in this module rather
than the others?", whereas a universe of all measured genes asks "is this set
enriched among the genes of this module at all?".
}
\value{
\code{module_enrichment()} returns a \code{data.table} with one row per tested
module / gene-set pair: \code{Cluster}, \code{Term}, \code{Overlap},
\code{ModuleSize}, \code{SetSize}, \code{Universe}, \code{FoldEnrichment},
\code{PValue}, \code{Padj}, and \code{Genes} (the shared genes), ordered by
module and p-value.

\code{enrichment_dotplot()} returns a ggplot object, also saved when
\code{file} is supplied.
}
\seealso{\code{\link{longitudinal_heatmap}}, \code{\link{multi_heatmap}}}
\examples{
data("heatmap_example", package = "ModuleViz")
obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
                            time_col = "TimeHr", group_col = "Condition",
                            k = 4, stability = FALSE)
genes <- rownames(heatmap_example$mat)
sets <- list(Early_response = genes[1:15], Mid_response = genes[16:30],
             Late_response = genes[31:45],
             Housekeeping = genes[seq(2, 60, by = 4)])
enr <- module_enrichment(obj, sets, min_size = 3)
head(enr)

enrichment_dotplot(enr)

## top terms beside each module of the heatmap
longitudinal_heatmap(obj, enrichment = enr)
}
