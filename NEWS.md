# ModuleViz 0.1.0

## New features

* `module_correlation_heatmap()` draws a gene-gene correlation heatmap for each
  temporal module, exposing modules that are averaging over more than one
  programme. Large modules are subsampled (`max_genes`, `select`); a `.pdf`
  destination collects every module as one multi-page file, any other extension
  writes one file per module.
* `module_correlation_summary()` returns the mean within-module correlation for
  every module without drawing anything.
* `longitudinal_heatmap()` and `dual_omics_heatmap()` gain `show_column_names`,
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
