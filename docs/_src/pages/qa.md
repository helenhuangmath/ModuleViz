## Installation

### `BiocManager::install("helenhuangmath/ModuleViz")` fails. What is missing?

Installing from a GitHub repository needs the `remotes` package. Run
`install.packages("remotes")` and try again. See
[Install ModuleViz](install.html).

### Can I install with `BiocManager::install("ModuleViz")`?

Not yet. `ModuleViz` is not on Bioconductor, so use the GitHub form.

## Input

### Do I need to z-score my matrix first?

No. `longitudinal_cluster()` z-scores each feature across samples. If the
matrix is already z-scored, pass `is_zscore = TRUE`.

### Which features should I cluster?

Select features before clustering, using differential analysis, a variance
filter, or another criterion that suits the experiment.
`select_variable_features()` keeps the most variable features of a matrix or a
`SummarizedExperiment`.

### My time variable is not numeric. Does that work?

Yes. `time_col` can be numeric, a factor (the level order is used), or
character (the order of first appearance is used).

### How are replicates handled?

Set `aggregate_replicates = TRUE` in `longitudinal_cluster()` to average
replicates within each group and time point before clustering.

## Clustering

### How do I choose the number of modules?

`choose_k()` reports the within-cluster sum of squares and the mean silhouette
width over a range of `k`, and suggests an elbow. Treat the elbow as a
diagnostic. Prefer the `k` that separates programmes you can interpret: a `k`
that merges two distinct modules is worse than one that splits a module in two.
Check the result with `cluster_stability()` and the
[module coherence](coherence.html) plots.

### Which clustering method should I use?

| `method` | Use it when |
|---|---|
| `"kmeans"` | The default. Fast. |
| `"hierarchical"` | You want a deterministic result. |
| `"pam"` | Outliers are a concern. |
| `"clara"` | There are tens of thousands of features. |
| `"fuzzy"` | You want a membership score per feature. |

Set `dist_method = "correlation"` when the shape of a trajectory matters more
than its amplitude.

### What do the module numbers mean?

Modules are ordered by the time point at which their centroid peaks. `Cluster1`
is the earliest-peaking module and the last cluster is the latest-peaking one.

### How do I read the stability values?

`cluster_stability()` subsamples features, reclusters, and reports the best
Jaccard overlap with the original module. High values mean the module survives
subsampling. Interpret values below about 0.6 with caution.

### A module looks like two patterns mixed together. What should I do?

Draw it with `module_correlation_heatmap()`. Visible sub-blocks suggest the
module should be split, which usually means increasing `k`.

## Multi-omics

### How are peaks matched to genes?

By row name. The first layer's row IDs are taken as gene symbols, and the other
layers use the text after the last underscore (`peak12_Pdcd1` maps to `Pdcd1`).
Supply `gene_of` for a different naming scheme.

### Several peaks map to one gene. Which one is shown?

`aggregate = "mean"` averages them, and `aggregate = "top"` keeps the strongest
feature.

### What happens to genes missing from one layer?

`gene_set = "intersect"` keeps only genes present in every layer.
`gene_set = "reference"` keeps every gene of the reference layer and draws the
missing rows in grey.

### Can modules be defined from all layers together?

Yes. `cluster_by = "joint"` clusters the concatenated per-gene profiles of all
layers, so a module groups genes whose expression and chromatin trajectories
agree. Otherwise `cluster_by` names the layer that defines the modules.

## Single-cell data

### Can I cluster single cells directly?

No. `ModuleViz` clusters features across an ordered set of samples, so cells
are first summarised with `pseudobulk_cells()`.

### My data has pseudotime, not time points.

Pass the pseudotime column as `time_col` and set `n_bins`. Cells are split into
equal-sized bins along it, and the bin index becomes the time axis.

### My data is already log-normalised.

Use `fun = "mean"` in `pseudobulk_cells()` to average values instead of summing
counts.

## Enrichment

### Which gene universe is used?

Supply `universe` (all measured genes) when the object holds only pre-selected
genes. Supply `gene_of` when the rows are peaks and not genes.

### Can I show results from another enrichment tool beside the heatmap?

Yes. Any data frame with `Cluster`, `Term`, and `Padj` columns can be passed to
the `enrichment` argument of `longitudinal_heatmap()`, `multi_heatmap()`, and
`dual_heatmap()`.

## Figures

### How do I change colours?

Use `zscore_palette` for the heatmap body (for example `"rdbu"` or
`"viridis"`), `cluster_palette` for the module bar, and `annotation_palette`
for the sample annotation. Use `annotation_colors` when exact metadata colours
are needed.

### How do I label specific genes?

Use `label_features` in `multi_heatmap()` and `dual_heatmap()`, or
`top_n_label` in `longitudinal_heatmap()` to label the strongest features of
each module. `resolve_labels()` shows which rows will be labelled.
