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

# Initial release

* Initial package release.
* Adds temporal module clustering, ordered longitudinal heatmaps, module
  trajectory plots, membership exports, dual-omics heatmaps, public real-data
  loading helpers, examples, and a tutorial vignette.
