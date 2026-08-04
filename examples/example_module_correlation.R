## =============================================================================
## Gene-gene correlation heatmaps, one per temporal module.
##
## longitudinal_cluster() tells you which features share a trajectory; it does
## not tell you how tightly they agree.  A module can look clean in the main
## heatmap and still be the average of two anti-correlated sub-programmes.
## Correlating every feature in a module against every other one exposes that:
## a coherent module is a solid block of high correlation, while a module with
## sub-structure breaks into visible blocks.
##
## Run from the package root:
##   Rscript examples/example_module_correlation.R
## =============================================================================

library(ModuleViz)

out_dir <- file.path("examples", "outputs", "module_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data("heatmap_example")

obj <- longitudinal_cluster(
  heatmap_example$mat,
  heatmap_example$meta,
  sample_col = "SampleID",
  time_col   = "TimeHr",
  group_col  = "Condition",
  k = 4,
  aggregate_replicates = TRUE,
  stability = TRUE,
  n_resample = 20
)
print(obj)

## ---------------------------------------------------------------------------
## 1. The one-number summary: how coherent is each module?
## ---------------------------------------------------------------------------
## Roughly:  > 0.7  tight, coherent module
##           0.3-0.7 loose; check the heatmap for sub-blocks
##           < 0.3  mostly noise, or two opposing programmes averaging out
cor_summary <- module_correlation_summary(obj)
print(cor_summary)

data.table::fwrite(cor_summary,
                   file.path(out_dir, "module_correlation_summary.txt"),
                   sep = "\t")

## ---------------------------------------------------------------------------
## 2. One correlation heatmap per module, collected into a single PDF
## ---------------------------------------------------------------------------
## A .pdf destination writes one page per module.  Any other extension writes
## one file per module, with the module name inserted before the extension.
cors <- module_correlation_heatmap(
  obj,
  max_genes    = 60,        # cap per module so the labels stay legible
  select       = "top",     # keep the strongest-signal features
  corr_palette = "rdbu",
  file   = file.path(out_dir, "module_gene_correlation.pdf"),
  title  = "Example time course",
  width  = 7,
  height = 6.5
)

for (nm in names(cors)) {
  cat(sprintf("%s: %d features (%d shown), mean r = %.3f\n",
              nm, cors[[nm]]$n_total, cors[[nm]]$n_shown, cors[[nm]]$mean_cor))
}

## ---------------------------------------------------------------------------
## 3. A single module to PNG
## ---------------------------------------------------------------------------
## The module name is appended to the filename, so this writes
## single_module_Cluster1.png rather than overwriting one file per module.
module_correlation_heatmap(
  obj,
  clusters  = "Cluster1",
  max_genes = 40,
  file   = file.path(out_dir, "single_module.png"),
  width  = 7.5,
  height = 7
)

## ---------------------------------------------------------------------------
## 4. Options worth knowing
## ---------------------------------------------------------------------------
## Spearman instead of Pearson, for monotone but non-linear agreement:
print(module_correlation_summary(obj, method = "spearman"))

## An unbiased sample of a large module rather than its strongest features:
invisible(module_correlation_heatmap(
  obj, clusters = "Cluster2", select = "random", max_genes = 30, seed = 42,
  file = file.path(out_dir, "random_sample.png"), width = 7, height = 6.5))

## Viridis instead of the diverging default, and a narrower colour range so
## differences among already-high correlations are visible:
invisible(module_correlation_heatmap(
  obj, clusters = "Cluster2", corr_palette = "viridis", corr_limits = c(0, 1),
  max_genes = 40,
  file = file.path(out_dir, "viridis.png"), width = 7, height = 6.5))

## The raw correlation matrices come back too, for downstream use:
cm <- cors[["Cluster1"]]$cor
cat(sprintf("\nCluster1 correlation matrix: %d x %d\n", nrow(cm), ncol(cm)))
print(round(cm[1:4, 1:4], 3))

cat("\noutputs in:", normalizePath(out_dir), "\n")
print(list.files(out_dir))
