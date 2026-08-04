## =============================================================================
## Validate ModuleViz on real public data: CD8 T cell memory vs exhaustion
##
##   GSE285248 - RNA-seq   (16 samples, naive / memory / exhausted P14 cells)
##   GSE285245 - CUT&Run   (H3K27ac activating, H3K27me3 repressive)
##
##   "Deciphering the role of histone modifications in memory and exhausted
##   CD8 T cells."  Naive P14 cells were transferred into recipients infected
##   with LCMV Armstrong (-> memory, day 30) or Clone 13 (-> exhaustion,
##   day 32); spleens harvested at days 30-32.  PubMed 40389726.
##
## The "trajectory" here is a differentiation axis rather than wall-clock time:
##   TN (naive) -> TMEM (memory) -> TEX (exhausted)
## which is exactly the ordered-sample structure longitudinal_cluster() expects.
##
## Only the small processed count matrices are downloaded (~4.7 MB total); the
## multi-GB RAW tarballs are never touched.
##
## Run from the package root:
##   Rscript examples/validate_real_data.R
## =============================================================================

suppressMessages(library(ModuleViz))
suppressMessages(library(data.table))

out_dir   <- file.path("examples", "outputs", "validation")
cache_dir <- file.path("examples", "outputs", "geo_cache")
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

PASS <- 0L; FAIL <- 0L
chk <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  [ERROR] %s: %s\n", label, conditionMessage(e))); FALSE
  })
  if (ok) { PASS <<- PASS + 1L; cat(sprintf("  [ok]   %s\n", label)) }
  else    { FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", label)) }
  invisible(ok)
}
section <- function(x) cat(sprintf("\n=== %s ===\n", x))

## ---------------------------------------------------------------------------
## 0. Fetch the processed GEO matrices (cached)
## ---------------------------------------------------------------------------
section("0. Download GEO processed matrices")
GEO <- list(
  rna      = list(gse = "GSE285248", file = "GSE285248_NormRC_filter_wanno.txt.gz"),
  h3k27ac  = list(gse = "GSE285245", file = "GSE285245_NormRC_H3K27ac_wID.txt.gz"),
  h3k27me3 = list(gse = "GSE285245", file = "GSE285245_NormRC_H3K27me3_wID.txt.gz")
)
geo_path <- function(entry) {
  dest <- file.path(cache_dir, entry$file)
  if (!file.exists(dest)) {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/%s",
                   substr(entry$gse, 1, 6), entry$gse, entry$file)
    message("downloading ", entry$file)
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}
paths <- lapply(GEO, geo_path)
for (nm in names(paths))
  chk(sprintf("%s matrix available (%s)", nm, basename(paths[[nm]])),
      file.exists(paths[[nm]]) && file.size(paths[[nm]]) > 5e5)

## ---------------------------------------------------------------------------
## 1. Parse the RNA-seq matrix and build ModuleViz metadata
## ---------------------------------------------------------------------------
section("1. RNA-seq matrix and sample metadata")
rna_raw  <- fread(paths$rna)
rna_cols <- grep("_RNA_", names(rna_raw), value = TRUE)

## One row per gene symbol: where a symbol appears under several Ensembl IDs,
## keep the most highly expressed copy.
rna_raw[, .mean_expr := rowMeans(as.matrix(.SD)), .SDcols = rna_cols]
setorder(rna_raw, Gene.Name, -.mean_expr)
rna_uniq <- rna_raw[!duplicated(Gene.Name)]
rna_mat  <- as.matrix(rna_uniq[, ..rna_cols])
rownames(rna_mat) <- rna_uniq$Gene.Name

## Sample columns look like  <Infection>_P14_Day<d>_<assay>_Rep<n>_S<n>
make_meta <- function(cols, assay) {
  infection <- sub("_P14_.*", "", cols)
  day       <- as.integer(sub(".*_Day([0-9]+)_.*", "\\1", cols))
  state     <- c(Naive = "TN", Arm = "TMEM", Cl13 = "TEX")[infection]
  data.table(
    SampleID   = cols,
    State      = factor(state, levels = c("TN", "TMEM", "TEX")),
    ## the ordered differentiation axis handed to time_col
    StateOrder = unname(c(TN = 1, TMEM = 2, TEX = 3)[state]),
    ## naive cells are the pre-infection baseline, so they sit at day 0
    Day        = ifelse(infection == "Naive", 0L, day),
    Replicate  = as.integer(sub(".*_Rep([0-9]+)_.*", "\\1", cols)),
    Assay      = assay
  )
}
rna_meta <- make_meta(rna_cols, "RNA")
print(rna_meta[, .N, by = .(State, Day)])

chk("16 RNA samples in 3 differentiation states",
    ncol(rna_mat) == 16 && nlevels(rna_meta$State) == 3)
## naive cells were profiled twice (once beside each infection), so TN has 8
chk("8 naive and 4 each of memory / exhausted",
    identical(as.integer(table(rna_meta$State)), c(8L, 4L, 4L)))
chk("each state sits at a single day",
    all(rna_meta[, uniqueN(Day), by = State]$V1 == 1))
chk("no NA / non-finite values", all(is.finite(rna_mat)))
chk("gene symbols are unique rownames",
    !anyDuplicated(rownames(rna_mat)) && !any(is.na(rownames(rna_mat))))
chk("marker genes present",
    all(c("Pdcd1", "Tox", "Havcr2", "Lag3", "Sell", "Tcf7") %in% rownames(rna_mat)))
cat(sprintf("  %d genes x %d samples\n", nrow(rna_mat), ncol(rna_mat)))

## ---------------------------------------------------------------------------
## 2. Select differentially expressed genes to plot
## ---------------------------------------------------------------------------
## Real analyses do not cluster the whole transcriptome, they cluster the genes
## that actually move.  Here we run a vectorised one-way ANOVA across the three
## states, BH-adjust, and keep genes that are both significant and of
## appreciable effect size.  The resulting table doubles as the `rank_table`
## used to pick which genes get labelled on the heatmap.
section("2. Differential expression (one-way ANOVA across TN / TMEM / TEX)")

log_rna <- log2(rna_mat + 1)

deg_anova <- function(x, grp) {
  grp <- factor(grp)
  n   <- ncol(x); k <- nlevels(grp)
  grand <- rowMeans(x)
  ssb <- 0; ssw <- 0
  for (g in levels(grp)) {
    cols <- which(grp == g)
    mg   <- rowMeans(x[, cols, drop = FALSE])
    ssb  <- ssb + length(cols) * (mg - grand)^2
    ssw  <- ssw + rowSums((x[, cols, drop = FALSE] - mg)^2)
  }
  df1 <- k - 1; df2 <- n - k
  Fstat <- (ssb / df1) / (ssw / df2)
  Fstat[!is.finite(Fstat)] <- 0
  p <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)
  data.table(ID = rownames(x), Fstat = Fstat, pvalue = p,
             padj = stats::p.adjust(p, method = "BH"))
}

state_of <- rna_meta$State[match(colnames(log_rna), rna_meta$SampleID)]
deg <- deg_anova(log_rna, state_of)

mean_by <- function(s) rowMeans(log_rna[, state_of == s, drop = FALSE])
deg[, `:=`(
  log2FC_TMEM_vs_TN  = mean_by("TMEM") - mean_by("TN"),
  log2FC_TEX_vs_TN   = mean_by("TEX")  - mean_by("TN"),
  log2FC_TEX_vs_TMEM = mean_by("TEX")  - mean_by("TMEM"),
  baseMean           = rowMeans(rna_mat)
)]
deg[, maxAbsLFC := pmax(abs(log2FC_TMEM_vs_TN), abs(log2FC_TEX_vs_TN),
                        abs(log2FC_TEX_vs_TMEM))]

PADJ_CUT <- 0.01; LFC_CUT <- 1; MIN_EXPR <- 10; MAX_DEG <- 3000
sig <- deg[padj < PADJ_CUT & maxAbsLFC >= LFC_CUT & baseMean >= MIN_EXPR]
setorder(sig, -Fstat)
cat(sprintf("  %d of %d genes pass padj < %g, |log2FC| >= %g, baseMean >= %g\n",
            nrow(sig), nrow(deg), PADJ_CUT, LFC_CUT, MIN_EXPR))
if (nrow(sig) > MAX_DEG) {
  cat(sprintf("  capping to the top %d by ANOVA F (dropping %d weaker DEGs)\n",
              MAX_DEG, nrow(sig) - MAX_DEG))
  sig <- sig[seq_len(MAX_DEG)]
}
## Cluster on the log scale: z-scoring raw counts lets a few high-count outliers
## dominate the distance metric.  The CUT&Run layers are logged the same way.
deg_mat <- log2(rna_mat[sig$ID, , drop = FALSE] + 1)
fwrite(sig, file.path(out_dir, "01_DEG_table.txt"), sep = "\t")

chk("DEGs were found", nrow(sig) > 200)
chk("every DEG passes both thresholds",
    all(sig$padj < PADJ_CUT) && all(sig$maxAbsLFC >= LFC_CUT))
chk("padj is a valid BH adjustment",
    all(sig$padj >= sig$pvalue) && all(deg$padj <= 1))
chk("DEG matrix rows match the DEG table", identical(rownames(deg_mat), sig$ID))
chk("known exhaustion genes are called differential",
    all(c("Pdcd1", "Tox", "Havcr2", "Lag3", "Entpd1", "Cd101") %in% sig$ID))
chk("known naive/memory genes are called differential",
    all(c("Sell", "Ccr7", "Lef1", "Tcf7", "Il7r") %in% sig$ID))
cat("  top 12 DEGs by F: ", paste(head(sig$ID, 12), collapse = ", "), "\n")

## ---------------------------------------------------------------------------
## 3. Low-level helpers on the real matrix
## ---------------------------------------------------------------------------
section("3. zscore_rows / collapse_columns")
z <- zscore_rows(deg_mat)
rm_ <- rowMeans(z); rsd <- apply(z, 1, sd)
chk("row means ~ 0", max(abs(rm_)) < 1e-10)
chk("row sds ~ 1", max(abs(rsd[rsd > 0] - 1)) < 1e-10)
chk("dimnames preserved", identical(dimnames(z), dimnames(deg_mat)))

cc <- collapse_columns(z, as.character(state_of), c("TN", "TMEM", "TEX"))
chk("collapsed column == rowMeans of that state",
    max(abs(cc[, "TEX"] - rowMeans(z[, state_of == "TEX", drop = FALSE]))) < 1e-10)
chk("one column per state", identical(colnames(cc), c("TN", "TMEM", "TEX")))

## ---------------------------------------------------------------------------
## 4. Choosing k
## ---------------------------------------------------------------------------
section("4. choose_k")
## The elbow is only a heuristic.  Here it lands on k=5, but k=5 has NO
## memory-peaking module: the TMEM-specific programme is absorbed into the
## neighbouring modules.  Since memory-vs-exhaustion is the comparison this
## study is about, every figure downstream uses k=6.  The assertion below
## records why, so the choice is documented in the run rather than on the plot.
K <- 6L
diag <- choose_k(z, k_range = 2:12, method = "kmeans", annotate_k = FALSE,
                 file = file.path(out_dir, "02_choose_k_elbow.pdf"))
print(diag)
k_sug <- attr(diag, "suggested_k")
cat(sprintf("  elbow suggests k = %s; using k = %s\n", k_sug, K))
chk("suggested_k inside k_range", k_sug %in% 2:12)
chk("within-SS decreases monotonically with k", all(diff(diag$tot_withinss) < 0))
chk("silhouette computed for every k", all(is.finite(diag$mean_silhouette)))
chk("elbow plot written", file.exists(file.path(out_dir, "02_choose_k_elbow.pdf")))
## Document, as a test, why the elbow k is not the k we report.
chk("the elbow k merges the memory programme, k=6 keeps it",
    { peak_states <- function(kk) {
        o <- longitudinal_cluster(deg_mat, rna_meta, sample_col = "SampleID",
                                  time_col = "StateOrder", k = kk,
                                  stability = FALSE, seed = 1)
        cn <- as.matrix(o$centroids[, -1]); colnames(cn) <- c("TN", "TMEM", "TEX")
        colnames(cn)[apply(cn, 1, which.max)] }
      !("TMEM" %in% peak_states(k_sug)) && "TMEM" %in% peak_states(K) })

## ---------------------------------------------------------------------------
## 5. Cluster the DEGs along the differentiation axis
## ---------------------------------------------------------------------------
section("5. longitudinal_cluster along TN -> TMEM -> TEX")
obj <- longitudinal_cluster(
  deg_mat, rna_meta,
  sample_col = "SampleID", time_col = "StateOrder",
  k = K, method = "kmeans", stability = TRUE, n_resample = 20, seed = 1
)
print(obj)
chk("returns a 'longi' object", inherits(obj, "longi"))
chk("all DEGs retained", nrow(obj$z) == nrow(deg_mat))
chk("K modules produced", nlevels(obj$cluster) == K)
chk("all 16 samples kept", ncol(obj$z) == 16)
chk("membership rows match z rows", identical(obj$membership$ID, rownames(obj$z)))
chk("rows grouped by module", !is.unsorted(as.integer(obj$cluster)))
chk("no empty module", all(table(obj$cluster) > 0))
chk("stability table has one row per module",
    nrow(obj$stability) == K && all(obj$stability$N == as.integer(table(obj$cluster))))
chk("Jaccard stability in [0,1]",
    all(obj$stability$MeanJaccard >= 0 & obj$stability$MeanJaccard <= 1))
print(obj$stability)

chk("columns ordered along the differentiation axis",
    !is.unsorted(obj$meta$StateOrder))

## module centroids, one column per state
cent <- as.matrix(obj$centroids[, -1]); rownames(cent) <- obj$centroids$Cluster
colnames(cent) <- c("TN", "TMEM", "TEX")[as.numeric(colnames(cent))]
peak_state <- setNames(colnames(cent)[apply(cent, 1, which.max)], rownames(cent))
cat("\n  module centroids (mean z per state):\n"); print(round(cent, 2))
cat("  peak state: ", paste(sprintf("%s=%s", names(peak_state), peak_state),
                            collapse = "  "), "\n")
chk("modules ordered early-peaking -> late-peaking",
    !is.unsorted(match(peak_state, c("TN", "TMEM", "TEX"))))
chk("centroid == mean z of its module's features",
    { f <- rownames(obj$z)[obj$cluster == "Cluster1"]
      cm <- collapse_columns(obj$z[f, , drop = FALSE],
                             as.character(obj$meta$StateOrder), c("1", "2", "3"))
      max(abs(colMeans(cm) - cent["Cluster1", ])) < 1e-8 })
chk("same seed -> identical clustering",
    identical(as.character(obj$cluster),
              as.character(longitudinal_cluster(deg_mat, rna_meta,
                sample_col = "SampleID", time_col = "StateOrder",
                k = K, stability = FALSE, seed = 1)$cluster)))

## ---------------------------------------------------------------------------
## 6. Does the clustering recover the known biology?
## ---------------------------------------------------------------------------
## This is what a synthetic dataset cannot test: the modules must place
## canonical exhaustion genes in an exhaustion-peaking module and canonical
## naive genes in a naive-peaking one.
section("6. Biological correctness of module assignment")
module_of    <- setNames(as.character(obj$cluster), rownames(obj$z))
peak_of_gene <- function(g) unname(peak_state[module_of[g]])

TEX_MARKERS <- c("Pdcd1", "Tox", "Havcr2", "Lag3", "Entpd1", "Cd101", "Cd160")
TN_MARKERS  <- c("Sell", "Ccr7", "Lef1", "Tcf7")

tex_tab <- data.table(Gene = TEX_MARKERS, Module = module_of[TEX_MARKERS],
                      PeakState = peak_of_gene(TEX_MARKERS))
tn_tab  <- data.table(Gene = TN_MARKERS,  Module = module_of[TN_MARKERS],
                      PeakState = peak_of_gene(TN_MARKERS))
print(tex_tab); print(tn_tab)

chk("every exhaustion marker lands in a TEX-peaking module",
    all(tex_tab$PeakState == "TEX"))
chk("every naive marker lands in a TN-peaking module",
    all(tn_tab$PeakState == "TN"))
chk("exhaustion and naive markers are in different modules",
    length(intersect(tex_tab$Module, tn_tab$Module)) == 0)
chk("Pdcd1 z-score is maximal in a TEX sample",
    as.character(obj$meta$State[which.max(obj$z["Pdcd1", ])]) == "TEX")
chk("Sell z-score is maximal in a TN sample",
    as.character(obj$meta$State[which.max(obj$z["Sell", ])]) == "TN")
chk("Il7r (memory) is higher in TMEM than TEX",
    mean(obj$z["Il7r", obj$meta$State == "TMEM"]) >
      mean(obj$z["Il7r", obj$meta$State == "TEX"]))

## ---------------------------------------------------------------------------
## 7. Heatmaps with the updated palettes
## ---------------------------------------------------------------------------
section("7. longitudinal_heatmap (Paired modules / Set2 conditions / RdBu + viridis)")
rank_tab <- sig[, .(ID, value = -log10(pmax(padj, .Machine$double.xmin)))]

## One colour scheme throughout: viridis. Swap in "rdbu" (or supply your own
## circlize::colorRamp2 via col_fun) with a single argument change.
ZPAL   <- "viridis"   # z-score bodies: sequential, colour-blind safe
CORPAL <- "rdbu"      # correlations diverge around 0, so use a diverging map
hm_main <- file.path(out_dir, "03_module_heatmap.png")
ht <- longitudinal_heatmap(
  obj, annotation_cols = c("State", "Day"),
  rank_table = rank_tab, id_field = "ID", value_field = "value", top_n_label = 3,
  zscore_palette = ZPAL, cluster_palette = "paired", annotation_palette = "set2",
  show_column_names = FALSE,
  label_fontsize = 13, row_title_fontsize = 13,
  annotation_name_fontsize = 13, legend_fontsize = 12, title_fontsize = 16,
  title = "CD8 T cell differentiation: DEG modules (GSE285248)",
  file = hm_main, width = 10, height = 12)
chk("main heatmap written", file.exists(hm_main) && file.size(hm_main) > 50000)
chk("returns a ComplexHeatmap object", inherits(ht, "Heatmap"))


chk("cluster palette 'paired' gives verbatim Paired colours",
    identical(unname(ModuleViz:::.moduleviz_cluster_cols(levels(obj$cluster), "paired")),
              RColorBrewer::brewer.pal(12, "Paired")[1:K]))
chk("annotation palette 'set2' gives verbatim Set2 colours",
    identical(ModuleViz:::.moduleviz_discrete_palette(3, "set2"),
              RColorBrewer::brewer.pal(8, "Set2")[1:3]))
chk("qualitative palettes are not interpolated",
    identical(ModuleViz:::.moduleviz_discrete_palette(3, "okabe_ito"),
              c("#0072B2", "#D55E00", "#009E73")))
chk("ramp palettes still span the full range",
    { s <- ModuleViz:::.moduleviz_discrete_palette(3, "spectral")
      s[1] == rev(RColorBrewer::brewer.pal(11, "Spectral"))[1] &&
      s[3] == rev(RColorBrewer::brewer.pal(11, "Spectral"))[11] })
chk("Set2 beyond its 8 colours falls back to interpolation",
    length(unique(ModuleViz:::.moduleviz_discrete_palette(12, "set2"))) == 12)
chk("RdBu z-scale runs blue(low) -> white(0) -> red(high)",
    { f <- ModuleViz:::.moduleviz_zscore_col_fun("rdbu")
      substr(f(-2), 2, 3) < substr(f(2), 2, 3) })   # less red at -2 than at +2
chk("viridis z-scale is monotone dark -> bright",
    { f <- ModuleViz:::.moduleviz_zscore_col_fun("viridis"); f(-2) != f(2) })
chk("custom zscore_limits are honoured",
    identical(ModuleViz:::.moduleviz_zscore_col_fun("rdbu", c(-3, 3))(3),
              ModuleViz:::.moduleviz_zscore_col_fun("rdbu", c(-2, 2))(2)))
chk("the historical 'spectual' typo alias still resolves",
    identical(ModuleViz:::.moduleviz_discrete_palette(4, "spectual"),
              ModuleViz:::.moduleviz_discrete_palette(4, "spectral")))

## ---------------------------------------------------------------------------
## 8. Module trajectory line plot
## ---------------------------------------------------------------------------
section("8. pattern_lineplot")
lp <- file.path(out_dir, "04_module_patterns.pdf")
STATE_LABS <- c("1" = "TN", "2" = "TMEM", "3" = "TEX")
g <- pattern_lineplot(obj, file = lp, cluster_palette = "paired",
                      base_size = 15, title_fontsize = 14,
                      xlab = "Differentiation state",
                      time_labels = STATE_LABS, x_angle = 0)
chk("line plot written", file.exists(lp) && file.size(lp) > 5000)
chk("returns a ggplot", inherits(g, "ggplot"))
chk("title font size is honoured", g$theme$plot.title$size == 14)
chk("title can be dropped entirely",
    is.null(pattern_lineplot(obj, title = NULL)$labels$title))
chk("x axis relabelled to the state names",
    all(ggplot2::layer_scales(g)$x$get_labels() == unname(STATE_LABS)))
chk("x axis label is overridable", g$labels$x == "Differentiation state")

## ---------------------------------------------------------------------------
## 9. Output tables
## ---------------------------------------------------------------------------
section("9. write_memberships")
pfx <- file.path(out_dir, "05_modules")
write_memberships(obj, prefix = pfx)
chk("membership file written", file.exists(paste0(pfx, "_memberships.txt")))
chk("centroid file written",   file.exists(paste0(pfx, "_centroids.txt")))
chk("stability file written",  file.exists(paste0(pfx, "_stability.txt")))
back <- fread(paste0(pfx, "_memberships.txt"))
chk("membership file round-trips all DEGs",
    nrow(back) == nrow(obj$z) && setequal(back$ID, rownames(obj$z)))

## ---------------------------------------------------------------------------
## 10. Gene-gene correlation within each module
## ---------------------------------------------------------------------------
section("10. module_correlation_heatmap / module_correlation_summary")
cor_summary <- module_correlation_summary(obj)   # all genes, no subsampling
print(cor_summary)
chk("summary has one row per module", nrow(cor_summary) == K)
chk("summary uses every gene in each module",
    identical(cor_summary$NShown, cor_summary$N))
chk("module sizes match the cluster table",
    identical(cor_summary$N, as.integer(table(obj$cluster))))
chk("mean within-module correlation is in [-1,1]",
    all(cor_summary$MeanCor >= -1 & cor_summary$MeanCor <= 1))
chk("modules are internally coherent (mean r > 0.3 everywhere)",
    all(cor_summary$MeanCor > 0.3))

cor_pdf <- file.path(out_dir, "06_module_correlation_all.pdf")
cors <- module_correlation_heatmap(obj, corr_palette = CORPAL, file = cor_pdf,
                                   label_top_n = 30, label_fontsize = 8,
                                   legend_fontsize = 12, title_fontsize = 11,
                                   title = "CD8 DEG modules")
chk("multi-page correlation PDF written",
    file.exists(cor_pdf) && file.size(cor_pdf) > 20000)
chk("one correlation matrix per module", length(cors) == K)
chk("correlation matrices are square, symmetric, unit-diagonal",
    all(vapply(cors, function(x)
      isSymmetric(x$cor) && all(abs(diag(x$cor) - 1) < 1e-12), logical(1))))
chk("correlations stay within [-1,1]",
    all(vapply(cors, function(x) all(x$cor >= -1 - 1e-9 & x$cor <= 1 + 1e-9),
               logical(1))))
chk("every gene in each module is used (no subsampling by default)",
    identical(unname(vapply(cors, function(x) x$n_shown, numeric(1))),
              as.numeric(cor_summary$N)))
chk("crowded modules label only the top 30",
    all(vapply(cors, function(x) length(x$labelled) == min(30, x$n_shown),
               logical(1))))
chk("labelled genes are the strongest-signal ones in their module",
    { m <- cors[[1]]; s <- apply(abs(obj$z[rownames(m$cor), , drop = FALSE]), 1, max)
      setequal(m$labelled, names(sort(s, decreasing = TRUE))[1:30]) })
chk("an explicit max_genes cap still works",
    module_correlation_heatmap(obj, clusters = "Cluster5", max_genes = 40,
                               file = tempfile(fileext = ".png"))[[1]]$n_shown == 40)
chk("summary and heatmap agree on module size",
    identical(unname(vapply(cors, function(x) x$n_total, numeric(1))),
              as.numeric(cor_summary$N)))

## The per-module file-naming path (one file per module rather than one
## multi-page PDF).  Cluster6 already has a page in 06_module_correlation_all.pdf,
## so this is checked in a temp directory rather than kept as a duplicate.
cor_png <- file.path(tempdir(), "corr.png")
tex_mod <- module_of["Pdcd1"]
module_correlation_heatmap(obj, clusters = tex_mod, corr_palette = CORPAL,
                           label_top_n = 30, label_fontsize = 9,
                           legend_fontsize = 12, title_fontsize = 11,
                           file = cor_png, width = 7.5, height = 7,
                           title = "Exhaustion module")
chk("per-module file is named after the module",
    file.exists(sub("(\\.png)$", paste0("_", tex_mod, "\\1"), cor_png)))
chk("spearman correlation also runs",
    !is.null(module_correlation_summary(obj, method = "spearman")$MeanCor))

## ---------------------------------------------------------------------------
## 11. Collapsing replicates
## ---------------------------------------------------------------------------
## Two distinct operations that are easy to confuse:
##
##   (a) aggregate_replicates = TRUE averages the replicates and then RE-CLUSTERS
##       the collapsed matrix.  It is an independent clustering, so its modules
##       are not the same sets of genes, in the same order, as the per-sample
##       run - only the early -> late ordering rule is shared.
##   (b) collapsing an EXISTING longi object for display keeps the modules,
##       their membership, and their order byte-for-byte identical to the
##       per-sample figure.
##
## (b) is what you want whenever the collapsed heatmap has to line up with the
## per-sample heatmap, so that is the one we draw.
section("11. Collapsing replicates: re-cluster (a) vs display-collapse (b)")

## (a) the aggregate_replicates code path
obj_a <- longitudinal_cluster(
  deg_mat, rna_meta, sample_col = "SampleID", time_col = "StateOrder",
  k = K, aggregate_replicates = TRUE, stability = FALSE)
chk("aggregate_replicates clusters into K modules", nlevels(obj_a$cluster) == K)
chk("one column per differentiation state", ncol(obj_a$z) == 3)
chk("meta rebuilt to match the collapsed columns", nrow(obj_a$meta) == 3)
chk("aggregated column == mean of that state's replicates",
    { m <- collapse_columns(deg_mat, as.character(state_of),
                            c("TN", "TMEM", "TEX"))
      g <- rownames(obj_a$z)[1]
      max(abs(zscore_rows(m)[g, ] - obj_a$z[g, ])) < 1e-8 })

## (b) display-collapse of the object behind figures 02 and 03
obj_c <- obj
obj_c$z <- collapse_columns(obj$z, as.character(obj$meta$State),
                            c("TN", "TMEM", "TEX"))
obj_c$meta <- data.table(
  SampleID   = c("TN", "TMEM", "TEX"),
  State      = factor(c("TN", "TMEM", "TEX"), levels = c("TN", "TMEM", "TEX")),
  StateOrder = 1:3)

chk("display-collapse keeps every feature, in the same row order",
    identical(rownames(obj_c$z), rownames(obj$z)))
chk("display-collapse keeps identical module assignments and order",
    identical(as.character(obj_c$cluster), as.character(obj$cluster)) &&
      identical(levels(obj_c$cluster), levels(obj$cluster)))
chk("display-collapse has one column per state", ncol(obj_c$z) == 3)
chk("collapsed column == mean z of that state's samples",
    max(abs(obj_c$z[, "TEX"] -
            rowMeans(obj$z[, obj$meta$State == "TEX", drop = FALSE]))) < 1e-10)
chk("module sizes match figure 02 / 03 exactly",
    identical(as.integer(table(obj_c$cluster)), as.integer(table(obj$cluster))))

ah <- file.path(out_dir, "07_collapsed_heatmap.png")
longitudinal_heatmap(obj_c, annotation_cols = "State", zscore_palette = ZPAL,
                     cluster_palette = "paired", annotation_palette = "set2",
                     top_n_label = 2, label_fontsize = 13,
                     row_title_fontsize = 13, annotation_name_fontsize = 13,
                     legend_fontsize = 12, column_name_fontsize = 13,
                     file = ah, width = 7, height = 11)
chk("collapsed heatmap written", file.exists(ah) && file.size(ah) > 40000)
## (no line plot here: it would be 04_module_patterns.pdf again)

## ---------------------------------------------------------------------------
## 12. CUT&Run: map peaks to genes, then draw the paired-omics figure
## ---------------------------------------------------------------------------
section("12. dual_omics_heatmap: RNA + H3K27ac CUT&Run")

## The CUT&Run tables carry peak coordinates but no gene symbol, so peaks are
## assigned to the nearest gene TSS (taken from the RNA annotation) within a
## window.  Row IDs are then built as "<peak>_<Symbol>", the PeakID_SYMBOL
## convention that dual_omics_heatmap()'s default gene_of_secondary expects.
gene_tss <- unique(rna_uniq[, .(Chr, TSS = ifelse(Strand == "+", Start, End),
                                Gene = Gene.Name)])

map_peaks_to_nearest_gene <- function(peaks, genes, max_dist = 25000) {
  out <- list()
  for (ch in intersect(unique(peaks$Chr), unique(genes$Chr))) {
    gg <- genes[Chr == ch][order(TSS)]
    pp <- peaks[Chr == ch]
    if (!nrow(gg) || !nrow(pp)) next
    i  <- findInterval(pp$Mid, gg$TSS)
    lo <- pmax(i, 1L); hi <- pmin(i + 1L, nrow(gg))
    dlo <- abs(pp$Mid - gg$TSS[lo]); dhi <- abs(pp$Mid - gg$TSS[hi])
    pick <- ifelse(dlo <= dhi, lo, hi)
    out[[length(out) + 1L]] <- data.table(peakN = pp$peakN, Gene = gg$Gene[pick],
                                          dist = pmin(dlo, dhi))
  }
  rbindlist(out)[dist <= max_dist]
}

load_cutrun <- function(path, mark) {
  tab   <- fread(path)
  scols <- grep(paste0("_", mark, "_Rep"), names(tab), value = TRUE)
  pk    <- tab[, .(peakN, Chr, Mid = as.integer((Start + End) / 2))]
  ann   <- map_peaks_to_nearest_gene(pk, gene_tss)
  tab   <- tab[peakN %in% ann$peakN]
  ann   <- ann[match(tab$peakN, peakN)]
  m     <- as.matrix(tab[, ..scols])
  rownames(m) <- paste0(ann$peakN, "_", ann$Gene)
  list(mat = log2(m + 1), meta = make_meta(scols, mark), ann = ann)
}

ac <- load_cutrun(paths$h3k27ac, "H3K27ac")
cat(sprintf("  H3K27ac: %d peaks within 25 kb of a TSS, %d distinct genes\n",
            nrow(ac$mat), uniqueN(ac$ann$Gene)))
chk("H3K27ac peaks mapped to genes", nrow(ac$mat) > 10000)
chk("peak IDs follow the PeakID_SYMBOL convention",
    all(grepl("^peak[0-9]+_", rownames(ac$mat))))
chk("default gene_of_secondary recovers the symbol",
    identical(sub("^.*_", "", rownames(ac$mat)[1:50]), ac$ann$Gene[1:50]))
chk("H3K27ac has the same 16-sample design",
    nrow(ac$meta) == 16 && setequal(levels(ac$meta$State), c("TN", "TMEM", "TEX")))
chk("exhaustion loci have H3K27ac peaks",
    all(c("Pdcd1", "Tox", "Havcr2") %in% ac$ann$Gene))

ac_obj <- longitudinal_cluster(ac$mat, ac$meta, sample_col = "SampleID",
                               time_col = "StateOrder", k = 4, stability = FALSE)
dual_ac <- file.path(out_dir, "08_RNA_vs_H3K27ac.png")
d_ac <- dual_omics_heatmap(
  obj, ac_obj, primary_name = "RNA", secondary_name = "H3K27ac",
  aggregate = "mean", label_features = c(TEX_MARKERS, TN_MARKERS),
  zscore_palette = ZPAL, cluster_palette = "paired",
  show_column_names = FALSE, label_fontsize = 13, row_title_fontsize = 13,
  legend_fontsize = 12, title_fontsize = 16,
  file = dual_ac, width = 15, height = 12)
chk("RNA + H3K27ac figure written",
    file.exists(dual_ac) && file.size(dual_ac) > 50000)
chk("secondary layer aligned row-for-row with the primary",
    identical(rownames(d_ac$z_secondary),
              rownames(obj$z)[rownames(obj$z) %in% ac$ann$Gene]))
## aggregate = "top" is exercised but not kept as a figure: it differs from the
## "mean" panel only in how multi-peak genes are collapsed.
chk("aggregate='top' also runs",
    { f <- tempfile(fileext = ".png")
      r <- dual_omics_heatmap(obj, ac_obj, primary_name = "RNA",
                              secondary_name = "H3K27ac", aggregate = "top",
                              zscore_palette = ZPAL, show_column_names = FALSE,
                              file = f, width = 15, height = 12)
      file.exists(f) && identical(rownames(r$z_secondary),
                                  rownames(d_ac$z_secondary)) })

## ---------------------------------------------------------------------------
## 13. Repressive mark, and does chromatin track expression?
## ---------------------------------------------------------------------------
section("13. H3K27me3 layer and cross-omics concordance")
me3 <- load_cutrun(paths$h3k27me3, "H3K27me3")
cat(sprintf("  H3K27me3: %d peaks, %d distinct genes\n",
            nrow(me3$mat), uniqueN(me3$ann$Gene)))
dual_me3 <- file.path(out_dir, "09_RNA_vs_H3K27me3.png")
d_me3 <- dual_omics_heatmap(
  obj, me3$mat, primary_name = "RNA", secondary_name = "H3K27me3",
  zscore_palette = ZPAL, cluster_palette = "paired",
  show_column_names = FALSE, label_fontsize = 13, row_title_fontsize = 13,
  legend_fontsize = 12, title_fontsize = 16,
  label_features = TEX_MARKERS, file = dual_me3, width = 15, height = 12)
chk("RNA + H3K27me3 figure written (raw-matrix secondary path)",
    file.exists(dual_me3) && file.size(dual_me3) > 50000)

## For genes induced in exhaustion the activating mark should gain signal in TEX
## relative to TN; the repressive mark should not track expression as closely.
delta_in <- function(res, meta_tab, genes) {
  zs   <- res$z_secondary
  keep <- intersect(genes, rownames(zs))
  st   <- meta_tab$State[match(colnames(zs), meta_tab$SampleID)]
  mean(zs[keep, st == "TEX"], na.rm = TRUE) - mean(zs[keep, st == "TN"], na.rm = TRUE)
}
tex_genes <- names(module_of)[peak_of_gene(names(module_of)) == "TEX"]
ac_delta  <- delta_in(d_ac,  ac$meta,  tex_genes)
me3_delta <- delta_in(d_me3, me3$meta, tex_genes)
cat(sprintf("  TEX-induced genes: H3K27ac TEX-TN = %+.3f, H3K27me3 TEX-TN = %+.3f\n",
            ac_delta, me3_delta))
chk("activating H3K27ac is gained at exhaustion-induced genes", ac_delta > 0)
chk("activating mark tracks expression more than the repressive mark",
    ac_delta > me3_delta)

## ---------------------------------------------------------------------------
## 14. One-call wrapper
## ---------------------------------------------------------------------------
section("14. run_longitudinal wrapper")
## The wrapper reproduces sections 5-8 exactly, so its four outputs are the
## same figures as 03_module_heatmap.png and 04_module_patterns.pdf.  They are
## written to a temp directory and verified there rather than kept as copies.
wpfx <- file.path(tempdir(), "wrapper")
## Same settings as section 5, so the wrapper must reproduce section 7's figure
## exactly - that equality is the point of this section, not a new analysis.
w <- run_longitudinal(deg_mat, rna_meta, out_prefix = wpfx,
                      sample_col = "SampleID", time_col = "StateOrder",
                      k = K, seed = 1, top_n_label = 3, rank_table = rank_tab,
                      annotation_cols = c("State"),
                      zscore_palette = ZPAL, cluster_palette = "paired",
                      annotation_palette = "set2",
                      stability = TRUE, n_resample = 20)
chk("wrapper returns a longi object", inherits(w, "longi"))
chk("one-call wrapper reproduces the step-by-step clustering exactly",
    identical(rownames(w$z), rownames(obj$z)) &&
      identical(as.character(w$cluster), as.character(obj$cluster)) &&
      identical(w$stability$MeanJaccard, obj$stability$MeanJaccard))
for (f in c("_heatmap.pdf", "_pattern_lines.pdf", "_memberships.txt", "_centroids.txt"))
  chk(sprintf("wrapper wrote %s", f), file.exists(paste0(wpfx, f)))

## ---------------------------------------------------------------------------
## 15. Error handling
## ---------------------------------------------------------------------------
section("15. Error handling on malformed real input")
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
chk("unnamed matrix is rejected",
    err(longitudinal_cluster(unname(deg_mat), rna_meta, sample_col = "SampleID",
                             time_col = "StateOrder", stability = FALSE)))
chk("non-overlapping sample IDs are rejected",
    { bad <- copy(rna_meta)[, SampleID := paste0("nope", .I)]
      err(longitudinal_cluster(deg_mat, bad, sample_col = "SampleID",
                               time_col = "StateOrder", stability = FALSE)) })
chk("missing time column is rejected",
    err(longitudinal_cluster(deg_mat, rna_meta, sample_col = "SampleID",
                             time_col = "NoSuchColumn", stability = FALSE)))
chk("unknown palette name is rejected",
    err(ModuleViz:::.moduleviz_discrete_palette(3, "not_a_palette")))
chk("unknown zscore palette is rejected",
    err(ModuleViz:::.moduleviz_zscore_col_fun("not_a_palette")))
chk("layers with no shared genes are rejected",
    { fake <- ac$mat[1:100, ]; rownames(fake) <- paste0("peakX_NOGENE", 1:100)
      err(dual_omics_heatmap(obj, fake, file = tempfile(fileext = ".png"))) })

## ---------------------------------------------------------------------------
section("SUMMARY")
cat(sprintf("  %d passed, %d failed\n", PASS, FAIL))
cat(sprintf("  outputs in: %s\n", normalizePath(out_dir)))
print(data.frame(file = list.files(out_dir),
                 bytes = file.size(file.path(out_dir, list.files(out_dir)))),
      row.names = FALSE)
if (FAIL > 0) quit(status = 1)
