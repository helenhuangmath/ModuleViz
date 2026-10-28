## =============================================================================
## Example: several data sets side by side (RNA + ATAC + H3K27ac + H3K27me3)
##
## Every layer is z-scored within itself, collapsed to one row per gene, and
## drawn as its own panel.  All panels share the module split and row order, so
## each row is one gene across every assay.  The modules come either from one
## reference layer (usually RNA) or from all layers clustered jointly.
## =============================================================================

library(data.table)

## Preferred after installing the package:
## install.packages("remotes")
## remotes::install_github("helenhuangmath/ModuleViz")
library(ModuleViz)

## During local development, you can instead source the R files:
# source("../R/longitudinal_heatmap.R")
set.seed(1)

out_dir <- "path/to/Results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- load the layers --------------------------------------------------------
## RNA: one row per gene (rownames = gene symbol), columns = samples.
## ATAC / CUT&RUN: rows = peaks named "PeakID_SYMBOL", columns = samples.
read_layer <- function(path) {
  tab <- fread(path)
  m   <- as.matrix(tab[, -1]); rownames(m) <- tab[[1]]
  m
}
rna      <- read_layer("path/to/RNA_NormRC_filtered.txt")
atac     <- read_layer("path/to/ATAC_NormRC_filtered.txt")
h3k27ac  <- read_layer("path/to/H3K27ac_NormRC_filtered.txt")
h3k27me3 <- read_layer("path/to/H3K27me3_NormRC_filtered.txt")

## Sample metadata.  One table covering every layer's sample IDs works; when the
## assays have their own sample sheets, pass a named list instead, e.g.
## meta = list(RNA = rna_meta, ATAC = atac_meta, H3K27ac = ac_meta, ...).
meta <- fread("path/to/sample_conditions.txt")
meta[, Genotype := sub("_.*$", "", Condition)]

layers <- list(RNA = rna, ATAC = atac, H3K27ac = h3k27ac, H3K27me3 = h3k27me3)

## ---- Option A: RNA defines the modules, chromatin follows -------------------
res <- multi_heatmap(
  layers,
  meta = meta, sample_col = "SampleID", time_col = "TimeHr",
  group_col = "Genotype",
  cluster_by = "RNA",
  k = 10, stability = TRUE,              # forwarded to longitudinal_cluster()
  aggregate = "mean",                    # collapse peaks per gene (or "top")
  label_features = c("Sox2", "Pax6", "Gfap"),
  file = file.path(out_dir, "RNA_ATAC_H3K27ac_H3K27me3_heatmap.pdf"),
  height = 14
)
multi_lineplot(res, file = file.path(out_dir, "multi_omics_patterns.pdf"))
write_memberships(res, prefix = file.path(out_dir, "multi_omics"))

## ---- Option B: cluster all layers together ----------------------------------
## Modules group genes whose expression AND chromatin trajectories agree.
## Each layer is down-weighted by 1/sqrt(#samples) so no layer dominates
## (balance_layers = FALSE turns this off).
joint <- multi_cluster(
  layers,
  meta = meta, sample_col = "SampleID", time_col = "TimeHr",
  group_col = "Genotype",
  cluster_by = "joint", k = 10
)
multi_heatmap(joint, file = file.path(out_dir, "multi_omics_joint.pdf"))
multi_lineplot(joint, file = file.path(out_dir, "multi_omics_joint_patterns.pdf"))

## ---- Option C: keep every RNA gene, even without a peak in some layers ------
## Missing cells are drawn in grey (na_col).
# multi_heatmap(layers, meta = meta, sample_col = "SampleID",
#                     time_col = "TimeHr", group_col = "Genotype", k = 10,
#                     gene_set = "reference", na_col = "grey85",
#                     file = file.path(out_dir, "multi_omics_all_rna_genes.pdf"))
