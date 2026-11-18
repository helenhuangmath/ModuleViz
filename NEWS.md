# ModuleViz 1.0.0

## New features

* **Single-cell data.** `pseudobulk_cells()` summarises a
  `SingleCellExperiment`, a `Seurat` object, or a (sparse) count matrix with
  per-cell metadata into pseudobulk profiles by sample, cell type, condition,
  and time point, returning a `SummarizedExperiment` that
  `longitudinal_cluster()` accepts directly. Counts are summed and converted to
  log2 CPM (`fun = "sum"`), or log-normalised values are averaged
  (`fun = "mean"`). `n_bins` cuts a continuous ordering such as pseudotime into
  bins for data without discrete time points.
* **Bioconductor containers.** `longitudinal_cluster()`, `run_longitudinal()`,
  and the layers of `multi_cluster()` accept a `SummarizedExperiment`; sample
  metadata is taken from `colData` and the assay is chosen with `assay`.
* **Module enrichment.** `module_enrichment()` tests every module of a `longi`
  or `multi_longi` object for over-representation of user-supplied gene sets (a
  named list, a two-column data frame, or a GMT file) with a hypergeometric
  test. `longitudinal_heatmap()`, `multi_heatmap()`, and `dual_heatmap()` gain
  `enrichment`, which prints the top terms beside each module
  (`enrichment_top_n`, `enrichment_padj`, `enrichment_fontsize`).
  `enrichment_dotplot()` draws the same table as a dot plot.
* **More clustering algorithms.** Besides `"kmeans"` (still the default),
  `"hierarchical"`, and `"pam"`, `method` now accepts `"clara"` (PAM on
  sub-samples, for tens of thousands of features) and `"fuzzy"` (fuzzy c-means,
  which adds a `Membership` score to the membership table). `choose_k()`
  supports all of them.
* **Publication-style figure defaults.** Every plot now draws with less text
  and smaller, consistent type (9-11 pt), sized for a journal column rather
  than a slide. Heatmaps print only the time point under each column
  (`column_labels = "time"`; `"sample"` restores the full sample names), show
  the group as a single top bar, drop the redundant module legend, and leave
  the stability bar off unless `show_stability = TRUE`. `multi_heatmap()` gives
  every layer the same panel width (`layer_widths`) and no longer adds an
  overall title; `dual_heatmap()` gains `annotation_cols` /
  `annotation_colors`. Line and dot plots use one clean theme (no grid, no
  strip boxes), have no title by default, and keep x labels horizontal.
  `module_correlation_heatmap()` labels 15 features by default and omits the
  dendrogram. Font sizes, titles, labels, and figure sizes can still be set
  through the existing arguments.
* **Trend plots.** `pattern_lineplot()` and `multi_lineplot()` stack the modules
  in a single column by default, matching the heatmap's row order; `ncol` wraps
  them into a grid in every case (it was previously ignored when the object had
  a group column). `multi_lineplot()` distinguishes groups by line type and
  gains `split_groups = TRUE` for the earlier one-column-per-group layout.
* `time_col` may be a factor (level order) or character (order of first
  appearance) as well as numeric.
* `select_variable_features()` keeps the most variable features of a matrix or
  `SummarizedExperiment`, and `simulate_single_cell_example()` generates a small
  single-cell time course for examples and tests.
* Heatmaps can be written to `.svg`, `.tiff`, and `.jpeg` as well as `.pdf` and
  `.png`.

## Bug fixes

* `aggregate_replicates = TRUE` no longer turns a non-numeric time column into
  `NA`, no longer breaks on group names containing `|`, and now keeps every
  metadata column that is constant within a group-by-time cell instead of
  dropping all but the sample, group, and time columns.
* `pattern_lineplot()` showed integer codes on the x axis for non-numeric time
  points without a group column, and `NA` with one; both now show the labels in
  order.
* Sample IDs that are numbers or factors were used as column positions rather
  than names when reordering the matrix.
* `dist_method = "correlation"` failed when a feature was constant (undefined
  correlation); such features are now treated as uncorrelated.
* `choose_k()` failed for a single candidate `k` or when no `k` in `k_range`
  was usable; its plot title no longer says "k-means" for other methods.
* Clustering with `k` equal to one or to the number of features, a single time
  point, or a module that is empty at one time point no longer errors.
* `resolve_labels()` matched explicit features by substring before trying a
  whole-token match, so `"Il2"` also labelled `Il2ra` and `Il21`; unmatched
  features are now reported.
* `multi_lineplot()` failed when layers stored time in different types.
* Random-number use no longer changes the session's RNG state: `seed` is
  applied locally, and `seed = NULL` defers to the caller's `set.seed()`.
* Features with missing or non-finite values are reported when dropped, and
  missing metadata columns, duplicated sample IDs, missing time values, and
  non-numeric input now produce clear errors.
* `.Rbuildignore` patterns were double-escaped and did not match, so `docs/`
  and local build artifacts were included in the source tarball.

## Package

* Meets Bioconductor submission checks: every exported function has a runnable
  example, the source tarball is below 5 MB, and dependencies on
  `SummarizedExperiment`, `S4Vectors`, `Matrix`, and `withr` are declared.
* `heatmap_example` is no longer lazy-loaded; use `data("heatmap_example")`.
* The deprecated alias `dual_omics_heatmap()` has been removed; use
  `dual_heatmap()`.

# ModuleViz 0.1.0

## New features

* `multi_cluster()`, `multi_heatmap()`, and `multi_lineplot()`
  extend the paired heatmap to any number of layers, for example
  RNA + ATAC + H3K27ac + H3K27me3. Layers can be `longi` objects or raw
  matrices with per-layer or shared metadata. Each layer is collapsed to one row
  per gene and drawn as its own panel with a shared module split and one z-score
  legend. `cluster_by` chooses which layer defines the modules, or `"joint"`
  clusters the concatenated profiles of all layers. `gene_set = "reference"`
  keeps genes missing from some layers and draws those cells in grey.
  `multi_lineplot()` draws each module's trajectory in every layer, with
  error bars across replicates.
* `dual_omics_heatmap()` is renamed `dual_heatmap()` and is now the two-layer
  case of `multi_heatmap()`. Its arguments and return value are unchanged, and
  the two layers now share a single z-score legend. The old name still works
  but is deprecated.
* Heatmaps use larger default text and are slightly narrower: gene labels are
  15 pt (from 10), module names 16 pt (from 10), sample names 13 pt,
  annotation names and legends 14 pt, and titles 17 pt. Default widths are 8 in
  for `longitudinal_heatmap()` (from 9), 10.5 in for `dual_heatmap()`
  (from 12), and 6.5 in for `module_correlation_heatmap()` (from 7), whose
  labels, legend, and title also grow to 12, 14, and 16 pt.
* Heatmap row titles now give only the module size (`n=`). Module stability
  remains available as the `Stability` bar in `longitudinal_heatmap()` and in
  the object's `stability` table.
* `multi_heatmap()` separates the data-type panels with a 10 mm gap,
  adjustable with `layer_gap`.
* `module_correlation_heatmap()` draws a square heatmap body by default, sized
  to fit the figure after the labels, title, and legend (`square = FALSE` lets
  it stretch as before).
* `write_memberships()` also accepts `multi_longi` objects.
* `module_correlation_heatmap()` draws a gene-gene correlation heatmap for each
  temporal module, exposing modules that are averaging over more than one
  programme. Large modules are subsampled (`max_genes`, `select`); a `.pdf`
  destination collects every module as one multi-page file, any other extension
  writes one file per module.
* `module_correlation_summary()` returns the mean within-module correlation for
  every module without drawing anything.
* `longitudinal_heatmap()` and `dual_heatmap()` gain `show_column_names`,
  for datasets whose sample IDs are long or uninformative and whose annotation
  bars already identify each column.
* All heatmaps gain `legend_fontsize` and `title_fontsize`, so legend and title
  text can be scaled up for print alongside the existing `label_fontsize`,
  `row_title_fontsize`, `column_name_fontsize`, and `annotation_name_fontsize`.
* `choose_k()` gains `chosen_k`, which marks the module count actually used as a
  solid line beside the dashed elbow, keeping the diagnostic figure consistent
  with the rest of an analysis when the reported `k` differs from the heuristic,
  and `annotate_k = FALSE` for a bare elbow curve with no markers or subtitle.
* `pattern_lineplot()` gains `group_colors`, a named vector of explicit colours
  for the group lines, matching what `annotation_colors` already does for
  heatmap annotation bars.
* `longitudinal_heatmap()` now warns when an entry of `annotation_cols` is not
  present in the object's metadata instead of dropping it silently. This bites
  after `aggregate_replicates = TRUE`, which rebuilds the metadata with only the
  sample, group, and time columns.
* `pattern_lineplot()` gains `title`, `title_fontsize`, `xlab`, `time_labels`,
  `x_angle`, `line_width`, and `point_size`. `time_labels` relabels the x axis
  when the ordering column is an integer rank (`1, 2, 3` -> `TN, TMEM, TEX`),
  and the title no longer inflates with `base_size`.
* `module_correlation_heatmap()` now uses every feature in a module by default
  (`max_genes = Inf`) and labels the strongest `label_top_n = 30` via link lines
  when printing all names would be unreadable. Ordering is computed on the
  correlation matrix directly rather than by re-correlating correlation
  profiles, which was needlessly slow for modules of ~1000 features.

## Documentation

* New main vignette, `vignette("ModuleViz")`: a full user guide covering input
  requirements, the clustering and module-ordering rules, stability, palettes,
  every plot type, paired omics, and an FAQ, with evaluated figures throughout.
* New vignette `vignette("cd8-exhaustion-published-data")`: an end-to-end
  walkthrough on published GEO data (GSE285248 RNA-seq + GSE285245 CUT&RUN).
* New `examples/example_module_correlation.R`.

## Initial release

* Initial package release.
* Adds temporal module clustering, ordered longitudinal heatmaps, module
  trajectory plots, membership exports, dual-omics heatmaps, public real-data
  loading helpers, examples, and a tutorial vignette.
