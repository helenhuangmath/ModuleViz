`ModuleViz` reads plain R objects and Bioconductor containers, and writes
tab-separated tables and figure files.

## Input

### Matrix and metadata table

The basic input is two objects:

- a numeric matrix with features (genes, peaks, regions) in rows and samples in
  columns, with row names and column names set;
- a metadata table with one row per sample.

The metadata table needs:

| Column | Argument | Required | Content |
|---|---|---|---|
| Sample ID | `sample_col` | yes | Matches the matrix column names. |
| Ordering | `time_col` | yes | Numeric, a factor (level order is used), or character (order of first appearance). |
| Group | `group_col` | no | Condition, genotype, treatment, cell type, or another factor. |

Values are z-scored by row inside `longitudinal_cluster()`. Pass
`is_zscore = TRUE` if the matrix is already z-scored.

### SummarizedExperiment

A `SummarizedExperiment` can be passed in place of the matrix and metadata
table to `longitudinal_cluster()`, `run_longitudinal()`, and as a layer of
`multi_cluster()`. Sample metadata is taken from `colData`, and `assay` selects
the assay to use.

### Single-cell objects

`pseudobulk_cells()` accepts a `SingleCellExperiment`, a `Seurat` object, or a
count matrix (dense or sparse) with per-cell metadata. It returns a
`SummarizedExperiment` of pseudobulk profiles. See
[Single-cell data](singlecell.html).

### Multi-omics layers

`multi_cluster()` and `multi_heatmap()` take a named list of layers, such as
`list(RNA = ..., ATAC = ..., H3K27ac = ...)`. Each layer is a `longi` object or
a raw matrix. Raw matrices take their sample metadata from `meta`, which is
either one shared table or a named list with one table per layer.

Layers are matched by gene. By default the first layer's row IDs are used as
gene symbols, and the other layers use the text after the last underscore, so
`peak12_Pdcd1` maps to `Pdcd1`. Use `gene_of` to supply a different mapping.

### Gene sets

`module_enrichment()` accepts gene sets in three forms:

- a named list, one character vector of genes per term;
- a two-column data frame (term, gene);
- the path to a GMT file.

## Objects

### `longi`

`longitudinal_cluster()` returns a `longi` object, a list holding the ordered
z-score matrix, module labels, the membership table, module centroids,
stability estimates, the reordered sample metadata, and the parameters of the
call. The elements used most often are:

| Element | Content |
|---|---|
| `obj$z` | z-scored matrix in heatmap row and column order |
| `obj$cluster` | module label of each feature |
| `obj$membership` | feature-to-module table |
| `obj$centroids` | mean profile of each module |
| `obj$stability` | per-module stability, when `stability = TRUE` |

Modules are named `Cluster1`, `Cluster2`, and so on, ordered by the time point
at which the module centroid peaks.

### `multi_longi`

`multi_cluster()` returns a `multi_longi` object: a list with `layers` (the
aligned z-score matrix and metadata of each layer), `cluster`, `genes`,
`membership`, `stability`, `reference` (the `longi` object the modules came
from), and `params`.

### Enrichment table

`module_enrichment()` returns a table with the columns `Cluster`, `Term`,
`Overlap`, `ModuleSize`, `SetSize`, `Universe`, `FoldEnrichment`, `PValue`,
`Padj`, and `Genes`. Any data frame with `Cluster`, `Term`, and `Padj` columns
can be passed to the `enrichment` argument of the heatmap functions, so results
from another enrichment tool can be shown the same way.

## Output files

### Module tables

`write_memberships(obj, prefix = "moduleviz")` writes tab-separated text files:

| File | Content |
|---|---|
| `moduleviz_memberships.txt` | one row per feature with its module |
| `moduleviz_centroids.txt` | mean profile of each module |
| `moduleviz_stability.txt` | per-module stability, when it was computed |

With `method = "fuzzy"`, the membership table also carries a `Membership` score
between 0 and 1.

### Figures

The heatmap and line-plot functions take a `file` argument, with `width` and
`height` to set the figure size. Heatmap functions draw to the current device when `file`
is not given. `pattern_lineplot()` and `multi_lineplot()` return `ggplot`
objects, which can also be saved with `ggplot2::ggsave()`.

`module_correlation_heatmap()` collects several modules as pages of one `.pdf`
file; with any other extension it writes one file per module.
