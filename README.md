# ModuleViz

`ModuleViz` is an R package for clustering and visualizing longitudinal
high-throughput omics data. It is designed for RNA-seq, ATAC-seq, CUT&Run,
pseudobulk single-cell summaries, and other matrix-like time-course assays.

Given a feature-by-sample matrix and a sample metadata table, `ModuleViz`
clusters features into temporal modules, orders those modules by trajectory, and
produces publication-ready heatmaps, module pattern plots, paired-omics
heatmaps, and within-module correlation diagnostics.

`ModuleViz` builds on
[`ComplexHeatmap`](https://bioconductor.org/packages/ComplexHeatmap/) for
heatmap rendering and uses standard R matrix/data-frame inputs so it can be
used in existing Bioconductor workflows.

## Features

| Task | Function |
|------|----------|
| Cluster genes, peaks, or other features into ordered temporal modules | `longitudinal_cluster()` |
| Draw clustered longitudinal heatmaps with sample annotations | `longitudinal_heatmap()` |
| Plot mean temporal trajectories for each module | `pattern_lineplot()` |
| Compare paired omics layers, such as RNA and ATAC/CUT&Run, aligned by gene | `dual_omics_heatmap()` |
| Inspect within-module gene-gene correlation structure | `module_correlation_heatmap()` |
| Summarize mean within-module correlations | `module_correlation_summary()` |
| Choose a cluster count with elbow and silhouette diagnostics | `choose_k()` |
| Estimate module stability by resampling | `cluster_stability()` |
| Export module memberships, centroids, and stability tables | `write_memberships()` |
| Run clustering, plotting, and export in one workflow | `run_longitudinal()` |

## Installation

Install the current development version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c("ComplexHeatmap", "circlize"))
install.packages(c("cluster", "data.table", "ggplot2", "RColorBrewer", "remotes"))

remotes::install_github("helenhuangmath/ModuleViz")
```

Optional packages used by vignettes and public-data examples:

```r
BiocManager::install(c("BiocStyle", "timecoursedata"))
install.packages(c("knitr", "rmarkdown", "testthat"))
```

For local development from this repository:

```r
devtools::load_all(".")
```

## Input Data

`ModuleViz` expects two inputs:

- A numeric matrix with features in rows and samples in columns. Rows may be
  genes, peaks, proteins, or any other measured feature. Feature IDs should be
  stored in `rownames(mat)`.
- A sample metadata table with one row per sample. The metadata must contain a
  sample ID column matching `colnames(mat)` and a numeric or ordered time column.
  A grouping column, such as condition, genotype, treatment, or cell state, is
  optional.

Values are z-scored by row automatically unless `is_zscore = TRUE` is supplied.
Feature selection should be done upstream; pass the package the features you
want to cluster and plot.

## Quick Start

```r
library(ModuleViz)

data("heatmap_example")
mat <- heatmap_example$mat
meta <- heatmap_example$meta

obj <- longitudinal_cluster(
    mat,
    meta,
    sample_col = "SampleID",
    time_col = "TimeHr",
    group_col = "Condition",
    k = 4,
    method = "kmeans",
    stability = TRUE
)

longitudinal_heatmap(
    obj,
    annotation_cols = c("Condition", "TimeHr"),
    top_n_label = 2,
    label_features = c("Gene001", "Gene016"),
    cluster_palette = "paired",
    file = "moduleviz_heatmap.pdf"
)

pattern_lineplot(
    obj,
    palette = "set2",
    file = "moduleviz_patterns.pdf"
)

write_memberships(obj, prefix = "moduleviz")
```

The same workflow can be run with one convenience function:

```r
run_longitudinal(
    mat,
    meta,
    out_prefix = "moduleviz",
    sample_col = "SampleID",
    time_col = "TimeHr",
    group_col = "Condition",
    k = 4,
    top_n_label = 2
)
```

## Example Figures

The bundled gallery script generates the main example outputs:

```r
source("examples/example_moduleviz_gallery.R")
```

### Longitudinal Module Heatmap

Features are clustered into temporal modules and ordered from early-peaking to
late-peaking patterns. Sample metadata can be displayed as annotation bars.

![Longitudinal clustered heatmap](inst/extdata/results/01_longitudinal_clustered_heatmap.png)

### Module Trajectory Plot

`pattern_lineplot()` summarizes the mean z-score trajectory for each module.

![Module line patterns](inst/extdata/results/02_module_line_patterns.png)

### Paired Omics Heatmap

`dual_omics_heatmap()` uses the primary layer, such as RNA, to define module
order and aligns a secondary layer, such as ATAC or CUT&Run, by gene.

![RNA ATAC side-by-side heatmap](inst/extdata/results/03_rna_atac_side_by_side_heatmap.png)

### Pseudobulk Single-Cell Example

The pseudobulk example shows the same workflow on grouped single-cell summaries.

```r
source("examples/example_single_cell_pseudobulk.R")
```

![Single-cell pseudobulk heatmap](inst/extdata/results/04_single_cell_pseudobulk_heatmap.png)

![Single-cell pseudobulk line patterns](inst/extdata/results/05_single_cell_pseudobulk_patterns.png)

### Published CD8 T Cell Exhaustion Example

The vignette `vignettes/cd8-exhaustion-published-data.Rmd` demonstrates an
end-to-end analysis using published RNA-seq and CUT&Run data.

![CD8 module heatmap](vignettes/figures/cd8_01_module_heatmap.png)

![CD8 module patterns](vignettes/figures/cd8_03_module_patterns.png)

![CD8 RNA versus H3K27ac heatmap](vignettes/figures/cd8_04_rna_vs_h3k27ac.png)

![CD8 module correlation heatmap](vignettes/figures/cd8_05_module_correlation.png)

## Choosing The Number Of Modules

`choose_k()` reports clustering diagnostics across candidate module counts and
stores the elbow estimate in the returned table.

```r
diag <- choose_k(
    zscore_rows(mat),
    k_range = 2:10,
    method = "kmeans",
    chosen_k = 4,
    file = "choose_k_elbow.pdf"
)

attr(diag, "suggested_k")
```

The elbow is a diagnostic. Choose a module count that is biologically
interpretable and stable for the dataset.

## Paired Omics Workflow

For paired assays, the primary layer drives clustering and row order. The
secondary layer is aligned by gene and collapsed to one row per gene when
multiple secondary features map to the same gene.

```r
rna_obj <- longitudinal_cluster(
    rna,
    meta,
    sample_col = "SampleID",
    time_col = "TimeHr",
    group_col = "Genotype",
    k = 10
)

dual_omics_heatmap(
    primary = rna_obj,
    secondary = atac,
    primary_name = "RNA",
    secondary_name = "ATAC",
    gene_of_secondary = function(id) sub("^.*_", "", id),
    aggregate = "mean",
    label_features = c("Sox2", "Pax6"),
    file = "rna_atac_heatmap.pdf"
)
```

See `examples/example_single_omics.R`, `examples/example_dual_omics.R`, and
`examples/example_module_correlation.R` for complete runnable templates.

## Public Time-Course Data

`load_real_timecourse_example()` loads supported public time-course datasets
from the optional Bioconductor `timecoursedata` package and returns a matrix and
metadata table ready for `ModuleViz`.

```r
available_real_timecourse_datasets()

example <- load_real_timecourse_example(
    dataset = "sorghum_leaf",
    genotype = "BT642",
    min_week = 3,
    top_n_features = 1000
)

obj <- longitudinal_cluster(
    example$mat,
    example$meta,
    sample_col = example$sample_col,
    time_col = example$time_col,
    group_col = example$group_col,
    k = 5,
    aggregate_replicates = TRUE,
    stability = TRUE
)
```

See `examples/example_real_timecoursedata.R` and
`vignettes/longitudinal-heatmap-tutorial.Rmd` for a full workflow.

## Vignettes

The package includes:

- `vignette("ModuleViz")`: complete user guide for clustering, heatmaps,
  palettes, paired omics, and module correlation plots.
- `vignette("cd8-exhaustion-published-data")`: published-data RNA/CUT&Run
  workflow.
- `vignette("longitudinal-heatmap-tutorial")`: tutorial using public
  time-course data.

## Package Checks

Before Bioconductor submission, run checks from a clean R/Bioconductor
environment with all required and suggested dependencies installed:

```bash
R CMD build .
R CMD check --no-manual ModuleViz_*.tar.gz
```

Run Bioconductor-specific checks:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("BiocCheck")
BiocCheck::BiocCheck(".")
```

## Citation

After installation, cite `ModuleViz` from R with:

```r
citation("ModuleViz")
```

## License

`ModuleViz` is distributed under the MIT license.
