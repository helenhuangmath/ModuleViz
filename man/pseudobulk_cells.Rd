\name{pseudobulk_cells}
\alias{pseudobulk_cells}
\title{Summarise single cells into pseudobulk profiles}
\description{
Group single cells by sample, cell type, condition, and time point (or by
pseudotime bin) and sum or average their expression. The result is a
\code{SummarizedExperiment} with one column per pseudobulk profile, ready for
\code{\link{longitudinal_cluster}} and the rest of the ModuleViz workflow.
}
\usage{
pseudobulk_cells(x, meta = NULL, time_col, sample_col = NULL,
  group_col = NULL, celltype_col = NULL, assay = NULL,
  fun = c("sum", "mean"), normalize = NULL, min_cells = 10,
  n_bins = NULL, features = NULL, prior_count = 1, scale_factor = 1e6)
}
\arguments{
\item{x}{A \code{SingleCellExperiment} or \code{SummarizedExperiment}, a
\code{Seurat} object, or a features-by-cells matrix (dense or sparse).}
\item{meta}{Per-cell metadata, one row per cell. Required for a matrix;
optional otherwise, where it replaces the object's own cell metadata.}
\item{time_col}{Metadata column giving the time point (or pseudotime) of each
cell.}
\item{sample_col}{Optional column identifying the biological sample or
replicate each cell came from.}
\item{group_col}{Optional column for a condition, genotype, or treatment.}
\item{celltype_col}{Optional column of cell-type labels. When given, every cell
type gets its own pseudobulk profiles.}
\item{assay}{Assay (or Seurat assay) to read. Defaults to \code{counts} for
\code{fun = "sum"} and \code{logcounts} for \code{fun = "mean"}.}
\item{fun}{\code{"sum"} (default; for raw counts) or \code{"mean"} (for values
that are already log-normalised).}
\item{normalize}{Add a log2 counts-per-million assay named \code{logcpm}.
Defaults to \code{TRUE} for \code{fun = "sum"}; ignored for
\code{fun = "mean"}.}
\item{min_cells}{Pseudobulk profiles built from fewer cells are dropped.}
\item{n_bins}{Optional number of equal-sized bins to cut a continuous
\code{time_col}, such as pseudotime, into. The bin index then becomes the time
and the mean of the original values is kept as \code{<time_col>_mean}.}
\item{features}{Optional subset of features to keep.}
\item{prior_count, scale_factor}{log2-CPM settings:
\code{log2(count / library_size * scale_factor + prior_count)}.}
}
\details{
Cells are grouped by every column supplied among \code{sample_col},
\code{celltype_col}, \code{group_col}, and \code{time_col}. Supplying
\code{sample_col} keeps biological replicates separate, which is what allows
\code{longitudinal_cluster()} to show replicate variation or average over it
with \code{aggregate_replicates = TRUE}.

ModuleViz clusters features across an ordered set of samples, so it does not
operate on individual cells: summarising first also removes most of the
drop-out noise of single-cell counts. For a Seurat object the \code{counts}
layer is read for \code{fun = "sum"} and the \code{data} layer for
\code{fun = "mean"}; this needs the SeuratObject package.
}
\value{
A \code{SummarizedExperiment}. With \code{fun = "sum"} it has the assays
\code{counts} and (unless \code{normalize = FALSE}) \code{logcpm}; with
\code{fun = "mean"} a single assay of averaged values. Its \code{colData} holds
the grouping columns, every other metadata column that is constant within a
profile, and \code{n_cells}.
}
\seealso{\code{\link{longitudinal_cluster}},
\code{\link{select_variable_features}},
\code{\link{simulate_single_cell_example}}}
\examples{
sce <- simulate_single_cell_example(cells_per_sample = 10)

## one profile per sample (replicates kept separate)
pb <- pseudobulk_cells(sce, time_col = "TimeHr", sample_col = "SampleID",
                       group_col = "Condition", min_cells = 5)
pb

obj <- longitudinal_cluster(pb, time_col = "TimeHr", group_col = "Condition",
                            k = 4, stability = FALSE)
obj

## one profile per cell type x sample
pb_ct <- pseudobulk_cells(sce, time_col = "TimeHr", sample_col = "SampleID",
                          group_col = "Condition", celltype_col = "CellType",
                          min_cells = 5)
dim(pb_ct)

## a continuous ordering (e.g. pseudotime) cut into bins
pb_bins <- pseudobulk_cells(sce, time_col = "TimeHr", n_bins = 4)
SummarizedExperiment::colData(pb_bins)
}
