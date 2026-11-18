# ModuleViz: temporal module discovery and visualization for longitudinal omics

<p class="badges">
<img alt="Version @VERSION@" src="https://img.shields.io/badge/version-@VERSION@-blue">
<img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green">
<img alt="R 4.3 or later" src="https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3">
<img alt="Built on ComplexHeatmap" src="https://img.shields.io/badge/built%20on-ComplexHeatmap-orange">
</p>

Latest release:

- GitHub: [helenhuangmath/ModuleViz](https://github.com/helenhuangmath/ModuleViz), version @VERSION@

## Introduction

`ModuleViz` is an R package for clustering and visualizing longitudinal
high-throughput omics data. It works on feature-by-sample matrices from
RNA-seq, ATAC-seq, CUT&RUN, pseudobulk single-cell summaries, and related
time-course or ordered-state experiments.

The question it answers is: which genes, peaks, or other features follow
similar patterns along an ordered biological axis? The axis can be clock time,
treatment dose, differentiation state, disease stage, exhaustion state, or any
other ordering you supply.

`ModuleViz` takes a numeric matrix and a sample metadata table, z-scores each
feature, clusters features into temporal modules, orders the modules from
early-peaking to late-peaking, and draws the result as heatmaps, module
trajectory plots, multi-omics heatmaps, and module coherence diagnostics.
Heatmaps are drawn with
[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/), so the
package fits into Bioconductor workflows.

![Longitudinal module heatmap: rows split into temporal modules, columns ordered by time](assets/results/01_longitudinal_clustered_heatmap.png)

## Changes for ModuleViz @VERSION@

- **Single-cell input.** `pseudobulk_cells()` summarises a
  `SingleCellExperiment`, a `Seurat` object, or a count matrix into pseudobulk
  profiles.
- **Bioconductor containers.** A `SummarizedExperiment` can be passed wherever
  a matrix and metadata table are expected.
- **Any number of omics layers.** `multi_cluster()`, `multi_heatmap()`, and
  `multi_lineplot()` align RNA, ATAC, histone marks, and other layers by gene.
- **Module enrichment.** `module_enrichment()` tests modules against gene sets
  and prints the top terms beside the heatmap.
- **More clustering algorithms.** `"clara"` and `"fuzzy"` join `"kmeans"`,
  `"hierarchical"`, and `"pam"`.

The full list is in the [release notes](changelog.html).

## Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

BiocManager::install("helenhuangmath/ModuleViz")
```

`ModuleViz` needs R 4.3 or later. See [Install ModuleViz](install.html) for
optional packages and a development setup.

## Usage

Cluster features into temporal modules and draw the heatmap:

```r
library(ModuleViz)
data("heatmap_example")

obj <- longitudinal_cluster(
    heatmap_example$mat, heatmap_example$meta,
    sample_col = "SampleID", time_col = "TimeHr", group_col = "Condition",
    k = 4
)
longitudinal_heatmap(obj, file = "moduleviz_heatmap.pdf")
```

Plot the mean trajectory of each module:

```r
pattern_lineplot(obj, file = "moduleviz_patterns.pdf")
```

Align several omics layers to the same modules:

```r
multi_heatmap(
    list(RNA = obj, ATAC = atac, H3K27ac = h3k27ac),
    meta = meta, sample_col = "SampleID", time_col = "TimeHr",
    file = "multi_heatmap.pdf"
)
```

Run clustering, both plots, and the table export in one call:

```r
run_longitudinal(mat, meta, out_prefix = "moduleviz",
                 sample_col = "SampleID", time_col = "TimeHr", k = 4)
```

The [Tutorial](tutorial.html) walks through the workflow step by step. The
exported functions are:

| Function | Description |
|---|---|
| [`longitudinal_cluster()`](functions.html#longitudinal_cluster) | Cluster features into ordered temporal modules. |
| [`longitudinal_heatmap()`](functions.html#longitudinal_heatmap) | Draw the annotated module heatmap. |
| [`pattern_lineplot()`](functions.html#pattern_lineplot) | Plot the mean trajectory of each module. |
| [`choose_k()`](functions.html#choose_k) | Elbow and silhouette diagnostics for the number of modules. |
| [`cluster_stability()`](functions.html#cluster_stability) | Estimate module reproducibility by resampling. |
| [`multi_cluster()`](reference/multi_cluster.html) | Align any number of omics layers by gene and cluster them. |
| [`multi_heatmap()`](reference/multi_heatmap.html) | Draw one heatmap panel per layer in a shared row order. |
| [`multi_lineplot()`](reference/multi_lineplot.html) | Plot each module's trajectory in every layer. |
| [`dual_heatmap()`](functions.html#dual_heatmap) | Two-layer shortcut, for example RNA beside ATAC. |
| [`module_enrichment()`](reference/module_enrichment.html) | Test modules for gene-set over-representation. |
| [`enrichment_dotplot()`](reference/module_enrichment.html) | Show enrichment results as a dot plot. |
| [`module_correlation_summary()`](functions.html#module_correlation_summary) | Summarise within-module correlation. |
| [`module_correlation_heatmap()`](functions.html#module_correlation_heatmap) | Draw feature-feature correlation within modules. |
| [`pseudobulk_cells()`](reference/pseudobulk_cells.html) | Summarise single cells into pseudobulk profiles. |
| [`select_variable_features()`](reference/select_variable_features.html) | Keep the most variable features before clustering. |
| [`write_memberships()`](functions.html#write_memberships) | Export module membership, centroid, and stability tables. |
| [`run_longitudinal()`](functions.html#run_longitudinal) | Run clustering, plotting, and export in one call. |
| [`zscore_rows()`](functions.html#zscore_rows) | Centre and scale each feature across samples. |
| [`resolve_labels()`](functions.html#resolve_labels) | Choose which rows are labelled on a heatmap. |
| [`collapse_columns()`](reference/longitudinal_cluster.html) | Average replicate samples within a group and time point. |
| [`load_real_timecourse_example()`](reference/load_real_timecourse_example.html) | Load a public time-course dataset ready for clustering. |
| [`available_real_timecourse_datasets()`](reference/load_real_timecourse_example.html) | List the supported public datasets. |
| [`simulate_single_cell_example()`](reference/simulate_single_cell_example.html) | Simulate a small single-cell time course for examples. |

## Documentation

- [Install ModuleViz](install.html)
- [Tutorial](tutorial.html): the core workflow from matrix to module tables
- [Function guide](functions.html): a worked example for each function
- [Multi-omics analysis](multiomics.html), [Module coherence](coherence.html),
  and [Single-cell data](singlecell.html)
- [Published data example](cd8.html): CD8 T cell exhaustion, RNA-seq with
  CUT&RUN
- [Input and output formats](formats.html)
- [Common Q & A](qa.html)
- [Vignettes](vignettes.html): the three rendered package vignettes
- [Function reference](reference/index.html): arguments and return values

## Contribute

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/helenhuangmath/ModuleViz/issues). A report
is easiest to act on when it includes a small reproducible example and the
output of `sessionInfo()`.

## Citation

After installation, get the citation from R:

```r
citation("ModuleViz")
```

> Huang H (2026). *ModuleViz: an R package for temporal module discovery and
> visualization in longitudinal omics data*. R package version @VERSION@,
> <https://github.com/helenhuangmath/ModuleViz>.

## Other useful links

- [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/): the
  heatmap engine behind every ModuleViz heatmap
- [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
  and [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/):
  accepted input containers
- [timecoursedata](https://bioconductor.org/packages/timecoursedata/): public
  time-course datasets used in the examples
- [Bioconductor](https://bioconductor.org/)
