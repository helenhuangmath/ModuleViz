# ModuleViz: an R package for temporal module discovery and visualization in longitudinal omics data

Hua Huang^1,*

^1 Department of Systems Pharmacology and Translational Therapeutics, University of Pennsylvania, Philadelphia, PA, USA

* Correspondence: helenhuang.math@gmail.com

## Abstract

### Summary

Longitudinal transcriptomic and epigenomic studies are increasingly used to
measure biological processes such as immune activation, differentiation,
disease progression, treatment response, and cellular exhaustion. These
experiments often produce feature-by-sample matrices with thousands of genes,
peaks, or other molecular features measured along an ordered biological axis.
Interpreting these data requires more than displaying a heatmap: analysts must
select dynamic features, cluster features with similar trajectories, order the
resulting modules, summarize module-level behavior, evaluate module coherence,
and export module assignments for downstream interpretation.

`ModuleViz` is an open-source R package that provides an integrated workflow
for temporal module discovery and visualization from longitudinal omics
matrices. The package z-scores features, clusters rows into temporal modules,
orders modules from earliest-peaking to latest-peaking patterns, generates
publication-ready heatmaps using `ComplexHeatmap`, summarizes module
trajectories with `ggplot2`, aligns paired omics layers such as RNA and
ATAC/CUT&RUN by gene, and exports module membership tables. Additional
diagnostics help users choose the number of modules and evaluate module
stability and within-module correlation structure. `ModuleViz` is designed for
bulk RNA-seq, pseudobulk single-cell RNA-seq, ATAC-seq, CUT&RUN, and related
quantitative omics data represented as standard R matrices.

### Availability and Implementation

`ModuleViz` is implemented in R and distributed under the MIT license. Source
code is available at <https://github.com/helenhuangmath/ModuleViz>. Rendered
documentation and tutorials are available at
<https://helenhuangmath.github.io/ModuleViz/articles/ModuleViz.html>.

### Keywords

longitudinal omics; time course; temporal modules; heatmap; RNA-seq; ATAC-seq;
CUT&RUN; Bioconductor; visualization

## Introduction

High-throughput omics experiments frequently measure ordered biological
processes rather than static end points. Examples include developmental time
courses, immune cell activation, dose-response experiments, disease-stage
comparisons, and differentiation or exhaustion trajectories. Such studies often
identify thousands of dynamic genes or regulatory elements. The biological
question is then not only which individual features change, but which groups of
features follow coordinated temporal programs.

General-purpose R visualization packages provide powerful building blocks for
matrix heatmaps and graphics. `ComplexHeatmap` offers highly flexible heatmap
construction for genomic matrices [@gu2016complexheatmap], while `ggplot2`
provides a grammar for statistical graphics [@wickham2016ggplot2]. The broader
R and Bioconductor ecosystems provide a reproducible foundation for genomic
software and high-throughput analysis [@r2026; @huber2015orchestrating].
However, longitudinal omics analysis typically requires a reproducible sequence
of operations that connects clustering, temporal ordering, module summaries,
paired-assay visualization, and exportable result tables. In practice, these
steps are often assembled through project-specific scripts, making analyses
harder to reuse, review, and reproduce.

`ModuleViz` [@moduleviz2026] addresses this gap by packaging a complete
matrix-based workflow for longitudinal omics visualization. The package accepts
a numeric matrix with features in rows and samples in columns, together with a
sample metadata table containing an ordered axis. It then identifies temporal
modules, orders them by their peak position along the axis, and returns both
graphics and structured tables for downstream analysis. The package
deliberately uses simple input data structures so it can be used after common
upstream workflows such as differential expression analysis, peak
quantification, pseudobulk aggregation, or user-defined feature filtering.

## Design and Implementation

`ModuleViz` is implemented as an R package with documented functions, testthat
unit tests, example datasets, vignettes, and rendered HTML tutorials. Its core
workflow is centered on a `longi` object returned by `longitudinal_cluster()`.
This object stores the z-scored matrix, ordered sample metadata, module labels,
row order, module centroids, stability estimates when requested, and a
membership table. Downstream plotting and export functions consume this object,
which keeps the workflow reproducible once clustering parameters have been
chosen.

### Inputs

The required inputs are deliberately minimal:

1. A numeric matrix with features in rows and samples in columns.
2. A metadata table with one row per sample.
3. A sample identifier column matching the matrix column names.
4. A numeric ordering column, such as time, dose, differentiation state, or
   disease stage.

An optional grouping column can represent condition, genotype, treatment,
infection, cell type, or another experimental factor. Feature selection is left
to the user or to upstream statistical tools, because the most appropriate
filter depends on the experiment. For example, users may cluster differentially
expressed genes, variable chromatin peaks, or a curated signature set.

### Temporal Module Discovery

The main analysis function, `longitudinal_cluster()`, performs row-wise
z-scoring unless the input has already been scaled. It supports k-means,
hierarchical clustering, and partitioning around medoids. Euclidean distance is
available for amplitude-aware clustering, while correlation distance can be used
when trajectory shape is more important than absolute magnitude.

After clustering, modules are ordered by their mean trajectory across the
ordered sample axis. For each module, `ModuleViz` calculates a centroid and
identifies the time point or ordered state at which the centroid reaches its
maximum. Modules are then renamed so that `Cluster1` corresponds to the
earliest-peaking module and later clusters correspond to later-peaking
programs. This ordering makes the resulting heatmap and module tables easier to
interpret than arbitrary cluster labels.

### Visualization

`longitudinal_heatmap()` draws the primary module heatmap using
`ComplexHeatmap`. Rows are split by temporal module, columns are ordered by the
sample axis, and user-selected metadata columns can be displayed as top
annotations. The function includes named palettes, z-score color maps,
adjustable figure dimensions, optional feature labels, and export to PDF or PNG.

`pattern_lineplot()` summarizes each module as a mean trajectory plot using
`ggplot2`. These plots provide a compact description of module behavior and are
useful for connecting heatmap patterns to biological interpretations such as
early induction, sustained activation, transient response, repression, or late
state-specific expression.

For paired-assay studies, `dual_omics_heatmap()` aligns a secondary assay to
primary modules by gene. For example, RNA-defined modules can be shown beside
ATAC-seq or CUT&RUN signal after mapping peaks to genes. The secondary layer
can be aggregated per gene using the mean or strongest matching feature. This
enables coordinated visualization of transcriptional and regulatory dynamics
without requiring users to manually synchronize row order across assays.

### Diagnostics and Export

`choose_k()` reports total within-cluster sum of squares and mean silhouette
width across candidate module numbers, providing an empirical diagnostic for
choosing `k`. `cluster_stability()` estimates the robustness of modules under
feature resampling using best Jaccard overlap with the original module labels.
`module_correlation_summary()` and `module_correlation_heatmap()` evaluate
within-module coherence by calculating feature-feature correlations within each
module. These diagnostics help users identify modules that are robust,
biologically interpretable, or potentially composed of multiple subprograms.

`write_memberships()` exports feature-to-module assignments, module centroids,
and stability summaries. These tables can be used as input to Gene Ontology,
pathway, motif, transcription factor, or other downstream enrichment analyses.
The convenience wrapper `run_longitudinal()` combines clustering, heatmap
generation, trajectory plotting, and table export for finalized workflows.

## Example Workflow

The package includes a bundled `heatmap_example` dataset and complete tutorials
under `vignettes/`. A typical analysis starts by loading a matrix and metadata
table, clustering features, and drawing the ordered heatmap:

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
    aggregate_replicates = TRUE,
    stability = TRUE
)

longitudinal_heatmap(
    obj,
    annotation_cols = c("Condition", "TimeHr"),
    top_n_label = 3,
    title = "Example longitudinal modules",
    file = "moduleviz_heatmap.pdf"
)

pattern_lineplot(obj, file = "moduleviz_patterns.pdf")
write_memberships(obj, prefix = "moduleviz")
```

The returned `longi` object contains the scaled matrix, module assignments,
ordered rows and columns, module centroids, stability estimates, and a
feature-level membership table. The same object can be passed directly to the
heatmap, line plot, paired-omics, correlation, and export functions.

## Case Study: CD8 T Cell Differentiation and Exhaustion

The repository includes a real-data vignette demonstrating `ModuleViz` on
published CD8 T cell differentiation data from the Gene Expression Omnibus. The
workflow uses RNA-seq data from GSE285248 to define temporal modules across an
ordered axis from naive cells to memory cells to exhausted cells, and aligns
H3K27ac CUT&RUN signal from GSE285245 to the RNA-defined modules by gene. This
case study illustrates an important design feature of `ModuleViz`: the ordered
axis does not need to be clock time. It can also represent an ordinal
biological state, such as differentiation stage, dose, disease grade, or
cell-state progression.

The vignette shows how module-level visualization can reduce a large gene-level
matrix to interpretable transcriptional programs. RNA modules reveal groups of
genes enriched at different positions along the naive-memory-exhausted axis,
and paired visualization enables inspection of whether chromatin activation
signals follow the same module order. The full workflow is available in
`vignettes/cd8-exhaustion-published-data.Rmd` and rendered at the project
website.

Suggested figures for a submitted manuscript are already present in the
repository:

- Figure 1: longitudinal clustered heatmap,
  `inst/extdata/results/01_longitudinal_clustered_heatmap.png`.
- Figure 2: module-average trajectory plot,
  `inst/extdata/results/02_module_line_patterns.png`.
- Figure 3: RNA/ATAC or RNA/CUT&RUN side-by-side heatmap,
  `inst/extdata/results/03_rna_atac_side_by_side_heatmap.png`.
- Figure 4: real-data CD8 module heatmap and paired H3K27ac visualization,
  `vignettes/figures/cd8_01_module_heatmap.png` and
  `vignettes/figures/cd8_04_rna_vs_h3k27ac.png`.

## Comparison With Existing Tools

`ModuleViz` is not intended to replace general-purpose heatmap or plotting
packages. Instead, it builds a higher-level longitudinal omics workflow around
them. Compared with manually combining `ComplexHeatmap`, clustering functions,
and custom plotting scripts, `ModuleViz` provides:

- a consistent object that stores clustered, ordered, and annotated results;
- temporal ordering of modules by centroid peak position;
- module-number and module-stability diagnostics;
- coordinated heatmap and trajectory outputs;
- paired-assay alignment by gene;
- module coherence summaries based on within-module correlation;
- exportable tables for downstream enrichment.

This scope makes `ModuleViz` complementary to the Bioconductor ecosystem. It
uses established R and Bioconductor graphics infrastructure while adding
domain-specific structure for ordered omics experiments.

## Reproducibility and Software Quality

The package contains documented exported functions, example data, vignettes,
rendered tutorials, and unit tests for core workflows. The test suite checks
clustering, `choose_k()` diagnostics, feature label resolution, and journal
palette options for trajectory plots. The package metadata declares runtime
dependencies in `Imports` and vignette or testing dependencies in `Suggests`.
Documentation is available both as R help pages and as rendered HTML tutorials.

For local validation, users can run:

```bash
R CMD build .
R CMD check --no-manual ModuleViz_*.tar.gz
```

and:

```r
testthat::test_dir("tests/testthat")
```

## Limitations and Future Work

`ModuleViz` is designed for matrix-based workflows and assumes that upstream
normalization and feature selection have already been performed. It does not
replace statistical differential expression, peak calling, motif enrichment, or
pathway enrichment tools. Instead, it provides the module discovery and
visualization layer that follows those analyses.

Future development could add direct support for common Bioconductor container
classes, interactive exploration, consensus clustering across datasets, richer
single-cell trajectory interfaces, and optional wrappers for downstream
functional enrichment. These additions would extend the current package while
preserving the simple matrix-and-metadata interface.

## Availability

Source code: <https://github.com/helenhuangmath/ModuleViz>

Documentation: <https://helenhuangmath.github.io/ModuleViz/articles/ModuleViz.html>

License: MIT

Package version described: 0.1.0

## Author Contributions

Hua Huang designed and implemented the software, prepared the documentation and
examples, and wrote the manuscript.

## Funding

Funding information should be added before journal submission.

## Conflict of Interest

The author declares no competing interests.

## Data Availability

The package includes a bundled example dataset for tutorials and tests. The
real-data vignette uses public Gene Expression Omnibus accessions GSE285248 and
GSE285245. Scripts and vignettes describing these analyses are included in the
repository.

## References

References are provided in `docs/paper/references.bib`.
