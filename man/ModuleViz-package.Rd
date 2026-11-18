\name{ModuleViz-package}
\alias{ModuleViz-package}
\alias{ModuleViz}
\title{Longitudinal Omics Heatmaps and Temporal Module Plots}
\description{
Cluster longitudinal omics features into temporal modules, draw ordered
ComplexHeatmap heatmaps, summarize module trajectories, test modules for gene-set
enrichment, and align paired omics layers such as RNA and ATAC by gene. Bulk
matrices, SummarizedExperiment objects, and single-cell data are supported.
}
\details{
Core workflow functions are \code{\link{longitudinal_cluster}},
\code{\link{longitudinal_heatmap}}, \code{\link{pattern_lineplot}}, and
\code{\link{write_memberships}}. Several omics layers are aligned by gene with
\code{\link{multi_cluster}} and \code{\link{multi_heatmap}}. Modules are tested
for gene-set enrichment with \code{\link{module_enrichment}}. Single-cell data
are summarised into pseudobulk profiles with \code{\link{pseudobulk_cells}}.
Public real-data examples can be loaded with
\code{\link{load_real_timecourse_example}}.
}
\value{
No return value; this page documents the package as a whole.
}
\author{Hua Huang}
\keyword{package}
