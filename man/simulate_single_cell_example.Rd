\name{simulate_single_cell_example}
\alias{simulate_single_cell_example}
\title{Simulate a small single-cell time course}
\description{
Generate single-cell counts whose sample-level averages follow the temporal
patterns of the bundled \code{\link{heatmap_example}} data: 60 genes, two
conditions, four time points, and two replicates, with several cell types per
sample. Intended for examples, tests, and tutorials.
}
\usage{
simulate_single_cell_example(cells_per_sample = 20,
  cell_types = c("T_cell", "Monocyte"), seed = 1)
}
\arguments{
\item{cells_per_sample}{Number of cells per sample and cell type.}
\item{cell_types}{Names of the simulated cell types.}
\item{seed}{RNG seed; \code{NULL} uses the current RNG state. The global
random-number stream is left unchanged.}
}
\value{
A \code{SingleCellExperiment} (a \code{SummarizedExperiment} if the
SingleCellExperiment package is not installed) with a \code{counts} assay and
the cell metadata columns \code{SampleID}, \code{Condition}, \code{TimeHr},
\code{Replicate}, and \code{CellType}.
}
\examples{
sce <- simulate_single_cell_example(cells_per_sample = 10)
sce
table(sce$CellType, sce$TimeHr)
}
