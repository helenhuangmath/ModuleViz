# ModuleViz

An installable R package for longitudinal / time-course omics analysis
(RNA-seq, ATAC-seq, Cut&Run, drug-treatment time courses, ...) built on
[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/).

Give it a matrix of features (genes or peaks) measured across an ordered set of
samples (time points / conditions) plus a small sample-metadata table, and it
will cluster the features into temporal modules and produce publication-ready
heatmaps and pattern plots.

## What it does

| # | Output | Function |
|---|--------|----------|
| 1 | Cluster features into temporal patterns/modules, with a choice of algorithm (k-means, hierarchical, PAM) and a **stability** assessment so modules are reproducible | `longitudinal_cluster()` |
| 2 | Ordered heatmap running from *high at the first time point* → *high at the last time point*, with metadata colour bars | `longitudinal_heatmap()` |
| 3 | Line/pattern plot of the mean trajectory per module | `pattern_lineplot()` |
| 4 | Write out cluster memberships (+ centroids + stability) | `write_memberships()` |
| 5 | Let the user pick specific features to label on the heatmap (plus auto top-N) | `label_features=` / `resolve_labels()` |
| 6 | Put **two data sets side by side** (RNA + ATAC, RNA + Cut&Run) so genes and peaks are shown together — a *separate* new figure | `dual_omics_heatmap()` |

Helpers: `choose_k()` (elbow + silhouette diagnostics to pick the number of
clusters), `cluster_stability()` (resampling Jaccard per module),
`run_longitudinal()` (one call that does cluster + heatmap + line plot + write).

The package code lives in [`R/`](R/); runnable templates are in
[`examples/`](examples/), including a real public-data example using
Bioconductor `timecoursedata`.

## Install

```r
install.packages(c("data.table", "circlize", "RColorBrewer", "ggplot2", "cluster"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "timecoursedata"))

install.packages("remotes")
remotes::install_github("helenhuangmath/Heatmap_Plot")
```

For local development from this folder:

```r
install.packages("devtools")
devtools::load_all(".")
```

## Package checks

Before submitting to GitHub, CRAN, or Bioconductor, run:

```bash
R CMD build .
env _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual ModuleViz_0.1.0.tar.gz
```

For Bioconductor-style checks:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BiocCheck")
BiocCheck::BiocCheck(".")
```

`timecoursedata` is optional and only needed for the public real-data examples.

## Inputs

- **Matrix** — features (rows) × samples (columns). Feature IDs as rownames
  (e.g. `"chr1:3046477-3046808_Xkr4"` or a gene symbol), or a data.table with an
  `id_col`. Values are z-scored by row automatically (`is_zscore = TRUE` to skip).
  *Feature selection (which genes/peaks are interesting) is done upstream* — pass
  in the matrix you want to plot.
- **Metadata** — one row per sample, with:
  - a sample-ID column matching the matrix columns (`sample_col`),
  - a **numeric time** column (`time_col`) used to order clusters and columns,
  - optionally a grouping column (`group_col`, e.g. genotype or treatment arm)
    that splits the heatmap columns and overlays lines in the pattern plot.

## Quick start

```r
library(ModuleViz)

data("heatmap_example")
mat <- heatmap_example$mat
meta <- heatmap_example$meta

# 1. cluster into temporal modules (ordered early-high -> late-high)
obj <- longitudinal_cluster(
  mat, meta,
  sample_col = "SampleID", time_col = "TimeHr", group_col = "Condition",
  k = 4, method = "kmeans",      # or "hierarchical" / "pam"
  stability = TRUE               # per-module reproducibility (mean Jaccard)
)
print(obj)

# 2. + 5. ordered heatmap with metadata bars and labelled genes
longitudinal_heatmap(
  obj,
  annotation_cols = c("Condition", "TimeHr"),
  top_n_label     = 2,                              # auto-label top genes...
  label_features  = c("Gene001", "Gene016"),        # ...plus chosen features
  file = "heatmap.pdf"
)

# 3. per-module pattern line plot
pattern_lineplot(obj, file = "pattern_lines.pdf")

# 4. cluster memberships (+ centroids + stability)
write_memberships(obj, prefix = "results")

# --- or everything at once ---
run_longitudinal(mat, meta, out_prefix = "results",
                 sample_col = "SampleID", time_col = "TimeHr",
                 group_col = "Condition", k = 4, top_n_label = 2)
```

## Example plot gallery

Run the bundled gallery script to generate examples of all three main outputs:

```r
source("examples/example_moduleviz_gallery.R")
```

It writes:

- `examples/outputs/moduleviz_gallery/01_longitudinal_clustered_heatmap.pdf`
- `examples/outputs/moduleviz_gallery/02_module_line_patterns.pdf`
- `examples/outputs/moduleviz_gallery/03_rna_atac_side_by_side_heatmap.pdf`
- `examples/outputs/moduleviz_gallery/moduleviz_example_memberships.txt`
- `examples/outputs/moduleviz_gallery/moduleviz_example_centroids.txt`
- `examples/outputs/moduleviz_gallery/moduleviz_example_stability.txt`

### Choosing the number of clusters

```r
diag <- choose_k(zscore_rows(mat), k_range = 2:20, method = "kmeans",
                 file = "elbow.pdf")
attr(diag, "suggested_k")   # distance-to-line elbow
```

### Two data sets side by side (requirement 6)

The primary layer (e.g. RNA) drives the clustering and row order; the secondary
layer (ATAC / Cut&Run) is aligned per gene and drawn alongside. Because one gene
can map to several peaks, the secondary layer is collapsed to one row per gene
(`aggregate = "mean"` or `"top"`).

```r
rna_obj <- longitudinal_cluster(rna, meta, sample_col = "SampleID",
                                time_col = "TimeHr", group_col = "Genotype", k = 10)

dual_omics_heatmap(
  primary = rna_obj, secondary = atac,
  primary_name = "RNA", secondary_name = "ATAC",
  gene_of_secondary = function(id) sub("^.*_", "", id),  # "PeakID_SYMBOL" -> "SYMBOL"
  aggregate = "mean",
  label_features = c("Sox2", "Pax6"),
  file = "RNA_ATAC_sidebyside.pdf"
)
```

See [`examples/example_single_omics.R`](examples/example_single_omics.R) and
[`examples/example_dual_omics.R`](examples/example_dual_omics.R) for full,
DESeq2-style templates.

## Real public-data example

The helper `load_real_timecourse_example()` loads public time-course data from
Bioconductor `timecoursedata` and returns a plain matrix + metadata table ready
for this package. The supported datasets are:

```r
available_real_timecourse_datasets()
```

Example using the Varoquaux 2019 sorghum leaf drought RNA-seq time course:

```r
library(ModuleViz)

example <- load_real_timecourse_example(
  dataset = "sorghum_leaf",
  genotype = "BT642",
  min_week = 3,
  top_n_features = 1000
)

diag <- choose_k(zscore_rows(example$mat), k_range = 2:10,
                 file = "choose_k_elbow.pdf")

obj <- longitudinal_cluster(
  example$mat,
  example$meta,
  sample_col = example$sample_col,
  time_col = example$time_col,
  group_col = example$group_col,
  k = attr(diag, "suggested_k"),
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 20
)

longitudinal_heatmap(obj, top_n_label = 3,
                     file = "sorghum_leaf_heatmap.pdf")
pattern_lineplot(obj, file = "sorghum_leaf_patterns.pdf")
write_memberships(obj, prefix = "sorghum_leaf")
```

See [`examples/example_real_timecoursedata.R`](examples/example_real_timecoursedata.R)
and the tutorial
[`vignettes/longitudinal-heatmap-tutorial.Rmd`](vignettes/longitudinal-heatmap-tutorial.Rmd).

Data source links:

- Bioconductor package page: <https://bioconductor.org/packages/timecoursedata/>
- Dataset documentation: <https://bioconductor.posit.co/packages/release/data/experiment/vignettes/timecoursedata/inst/doc/documentation.html>
- Sorghum tutorial context in `moanin`: <https://bioconductor.org/packages/release/bioc/vignettes/moanin/inst/doc/documentation.html>

## Notes on design

- **Temporal ordering.** Each module's centroid across the ordered time points is
  computed; modules are ranked by the time point where they peak (ties broken by
  overall slope), so the heatmap reads top-to-bottom from early-peaking to
  late-peaking. Rows within a module are ordered by strongest signal (max |z|).
- **Stability.** `cluster_stability()` sub-samples the features, re-clusters, and
  reports the mean maximum Jaccard overlap per module (à la `fpc::clusterboot`).
  Values near 1 mean a reproducible module; below ~0.6, interpret with caution.
  Stability is shown in the heatmap row titles (`J=…`) and as a left annotation.
- **Algorithms.** k-means (Hartigan–Wong, robust to empty clusters),
  hierarchical (`hclust`/`cutree`, default `ward.D2`), or PAM (`cluster::pam`).
  Distance can be Euclidean or `1 - correlation` (`dist_method = "correlation"`).

---

<details>
<summary>Other quick ways to plot a heatmap (reference snippets)</summary>

**pheatmap**
```r
pheatmap(scaled_mat, scale = "none", kmeans_k = 5, show_rownames = FALSE)
```

**ComplexHeatmap**
```r
Heatmap(scaled_mat, show_row_names = FALSE)
```

**heatmaply** — <https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html>
```r
heatmaply(as.matrix(dt[, 1:12]), scale = "row",
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"),
          file = "heatmap.pdf", height = 900, showticklabels = c(TRUE, FALSE))
```

**ggheatmap**
```r
ggheatmap(dt, scale = "column", row_side_colors = dt[, c("cyl", "gear")])
```
</details>
