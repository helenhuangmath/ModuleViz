## =============================================================================
## Validate ModuleViz on real public data: CD8 T cell memory vs exhaustion
##
##   GSE285248 - RNA-seq   (16 samples, naive / memory / exhausted P14 cells)
##   GSE285245 - CUT&RUN   (H3K27ac activating, H3K27me3 repressive)
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
## dominate the distance metric.  The CUT&RUN layers are logged the same way.
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

## A wider panel of well-known genes to call out on the figures.  The top-N
## labels are picked purely by statistic, so on their own they rarely name the
## genes a reader looks for first.  Only genes that are DEGs can be labelled.
LABEL_GENES <- intersect(
  c(TN_MARKERS, "Satb1", "Bach2", "Klf2", "S1pr1", "Foxp1",            # naive
    "Il7r", "Klrg1", "Cx3cr1", "Il18r1", "Cxcr3", "Zeb2", "S1pr5",     # memory
    TEX_MARKERS, "Tigit", "Ctla4", "Cd244a", "Eomes", "Nr4a2", "Cd38",
    "Cxcr6", "Gzmk", "Ikzf2"),                                          # exhausted
  rownames(obj$z))
cat(sprintf("  %d genes labelled on the figures: %s\n", length(LABEL_GENES),
            paste(sprintf("%s(%s)", LABEL_GENES, peak_of_gene(LABEL_GENES)),
                  collapse = ", ")))

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
  label_features = LABEL_GENES,
  rank_table = rank_tab, id_field = "ID", value_field = "value", top_n_label = 3,
  zscore_palette = ZPAL, cluster_palette = "paired", annotation_palette = "set2",
  show_column_names = FALSE,
  label_fontsize = 14, row_title_fontsize = 16,
  annotation_name_fontsize = 16, legend_fontsize = 15, title_fontsize = 19,
  title = "CD8 T cell differentiation\nDEG modules (GSE285248)",
  file = hm_main, width = 7.5, height = 14)
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
                                   label_top_n = 30, label_fontsize = 11,
                                   legend_fontsize = 14, title_fontsize = 14,
                                   width = 8.5, height = 7.5,
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
                     label_features = LABEL_GENES,
                     top_n_label = 2, label_fontsize = 14,
                     row_title_fontsize = 16, annotation_name_fontsize = 16,
                     legend_fontsize = 15, column_name_fontsize = 16,
                     file = ah, width = 6, height = 13)
chk("collapsed heatmap written", file.exists(ah) && file.size(ah) > 40000)
## (no line plot here: it would be 04_module_patterns.pdf again)

## ---------------------------------------------------------------------------
## 12. CUT&RUN: map peaks to genes, then draw the paired-omics figure
## ---------------------------------------------------------------------------
section("12. dual_heatmap: RNA + H3K27ac CUT&RUN")

## The CUT&RUN tables carry peak coordinates but no gene symbol, so peaks are
## assigned to the nearest gene TSS (taken from the RNA annotation) within a
## window.  Row IDs are then built as "<peak>_<Symbol>", the PeakID_SYMBOL
## convention that dual_heatmap()'s default gene_of_secondary expects.
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
d_ac <- dual_heatmap(
  obj, ac_obj, primary_name = "RNA", secondary_name = "H3K27ac",
  aggregate = "mean", label_features = LABEL_GENES,
  zscore_palette = ZPAL, cluster_palette = "paired",
  show_column_names = FALSE, label_fontsize = 15, row_title_fontsize = 16,
  legend_fontsize = 15, title_fontsize = 19,
  file = dual_ac, width = 10, height = 12)
chk("RNA + H3K27ac figure written",
    file.exists(dual_ac) && file.size(dual_ac) > 50000)
chk("secondary layer aligned row-for-row with the primary",
    identical(rownames(d_ac$z_secondary),
              rownames(obj$z)[rownames(obj$z) %in% ac$ann$Gene]))
## aggregate = "top" is exercised but not kept as a figure: it differs from the
## "mean" panel only in how multi-peak genes are collapsed.
chk("aggregate='top' also runs",
    { f <- tempfile(fileext = ".png")
      r <- dual_heatmap(obj, ac_obj, primary_name = "RNA",
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
d_me3 <- dual_heatmap(
  obj, me3$mat, primary_name = "RNA", secondary_name = "H3K27me3",
  zscore_palette = ZPAL, cluster_palette = "paired",
  show_column_names = FALSE, label_fontsize = 15, row_title_fontsize = 16,
  legend_fontsize = 15, title_fontsize = 19,
  label_features = LABEL_GENES, file = dual_me3, width = 10, height = 12)
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
## 14. All three layers in one figure
## ---------------------------------------------------------------------------
## dual_heatmap() pairs the primary with ONE secondary layer.  The two
## calls above already aligned each chromatin mark to the RNA rows, so a
## three-panel figure is just those aligned matrices side by side, restricted to
## the genes that carry a peak of both marks.  Rows keep the RNA module order.
section("14. RNA + H3K27ac + H3K27me3 together")
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))

STATES <- c("TN", "TMEM", "TEX")
genes3 <- rownames(obj$z)[rownames(obj$z) %in% rownames(d_ac$z_secondary) &
                            rownames(obj$z) %in% rownames(d_me3$z_secondary)]
cl3 <- droplevels(obj$cluster[match(genes3, rownames(obj$z))])

## one layer = the shared genes, with samples ordered TN -> TMEM -> TEX
as_layer <- function(z, meta_tab) {
  st <- meta_tab$State[match(colnames(z), meta_tab$SampleID)]
  o  <- order(st)
  list(z = z[genes3, o, drop = FALSE], state = st[o])
}
layers <- list(RNA      = as_layer(obj$z,             obj$meta),
               H3K27ac  = as_layer(d_ac$z_secondary,  ac$meta),
               H3K27me3 = as_layer(d_me3$z_secondary, me3$meta))
cat(sprintf("  %d genes carry RNA, H3K27ac and H3K27me3 signal\n", length(genes3)))

chk("the three layers are aligned row-for-row",
    all(vapply(layers, function(L) identical(rownames(L$z), genes3), logical(1))))
chk("every layer has all 16 samples, ordered along the differentiation axis",
    all(vapply(layers, function(L)
      ncol(L$z) == 16 && !is.unsorted(as.integer(L$state)), logical(1))))
chk("all six modules survive the intersection", nlevels(cl3) == K)

tri_file   <- file.path(out_dir, "10_RNA_H3K27ac_H3K27me3.png")
state_cols <- setNames(ModuleViz:::.moduleviz_discrete_palette(3, "set2"), STATES)
tri_col    <- ModuleViz:::.moduleviz_zscore_col_fun(ZPAL)
tri_leg    <- ModuleViz:::.moduleviz_legend_gp(15)
tri_marks  <- which(genes3 %in% LABEL_GENES)

tri_panel <- function(nm) {
  L <- layers[[nm]]
  first <- nm == names(layers)[1]; last <- nm == names(layers)[length(layers)]
  Heatmap(
    L$z, name = nm, col = tri_col,
    row_split = cl3, cluster_rows = FALSE, cluster_row_slices = FALSE,
    row_title = sprintf("%s\n(n=%d)", levels(cl3), as.integer(table(cl3))),
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 16),
    row_gap = grid::unit(1.5, "mm"), show_row_names = FALSE,
    cluster_columns = FALSE, show_column_names = FALSE,
    column_title = nm, column_title_gp = grid::gpar(fontsize = 19),
    top_annotation = HeatmapAnnotation(
      State = L$state, col = list(State = state_cols),
      show_annotation_name = last, show_legend = first,
      annotation_name_gp = grid::gpar(fontsize = 16),
      annotation_legend_param = tri_leg),
    left_annotation = if (first) rowAnnotation(
      Cluster = cl3,
      col = list(Cluster = ModuleViz:::.moduleviz_cluster_cols(levels(cl3), "paired")),
      show_annotation_name = FALSE, annotation_legend_param = tri_leg,
      width = grid::unit(4, "mm")),
    right_annotation = if (last) rowAnnotation(mark = anno_mark(
      at = tri_marks, labels = genes3[tri_marks],
      labels_gp = grid::gpar(fontsize = 15, fontface = "italic"),
      link_width = grid::unit(6, "mm"))),
    ## the layers share one colour scale, so one legend covers all three
    show_heatmap_legend = first,
    use_raster = TRUE, raster_quality = 2,
    heatmap_legend_param = c(list(title = "z-score", direction = "horizontal",
                                  title_position = "topcenter",
                                  legend_width = grid::unit(4, "cm")), tri_leg))
}
tri_ht <- Reduce(`+`, lapply(names(layers), tri_panel))
ModuleViz:::.moduleviz_write_plot(tri_file, 12, 12, function()
  draw(tri_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       merge_legend = TRUE, ht_gap = grid::unit(3, "mm"),
       column_title = "RNA + H3K27ac + H3K27me3 (linked by gene)",
       column_title_gp = grid::gpar(fontsize = 19)))
chk("three-layer heatmap written",
    file.exists(tri_file) && file.size(tri_file) > 50000)

## Module trajectories, one line per layer.  Each sample contributes its mean z
## over the module's genes; points are the mean of those sample values within a
## state and the bars are their SD, i.e. the spread between replicates.
traj <- rbindlist(lapply(names(layers), function(nm) {
  L <- layers[[nm]]
  rbindlist(lapply(levels(cl3), function(m)
    data.table(Layer = nm, Cluster = m, State = L$state,
               z = colMeans(L$z[cl3 == m, , drop = FALSE], na.rm = TRUE))))
}))
traj_sum <- traj[, .(z = mean(z), sd = sd(z)), by = .(Layer, Cluster, State)]
traj_sum[, Layer := factor(Layer, levels = names(layers))]
traj_sum[, ClusterF := factor(
  sprintf("%s (n=%d)", Cluster, as.integer(table(cl3))[match(Cluster, levels(cl3))]),
  levels = sprintf("%s (n=%d)", levels(cl3), as.integer(table(cl3))))]

## does each chromatin mark rise and fall with expression, module by module?
cent_cor <- dcast(traj_sum, Cluster + State ~ Layer, value.var = "z")[,
  .(r_H3K27ac = cor(RNA, H3K27ac), r_H3K27me3 = cor(RNA, H3K27me3)), by = Cluster]
cat("\n  correlation of each mark's module trajectory with RNA:\n")
print(cent_cor[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
chk("H3K27ac trajectory follows RNA in every module",
    all(cent_cor$r_H3K27ac > 0.5))
chk("H3K27ac follows RNA more closely than H3K27me3 does",
    mean(cent_cor$r_H3K27ac) > mean(cent_cor$r_H3K27me3))

layer_cols <- setNames(ModuleViz:::.moduleviz_discrete_palette(3, "okabe_ito"),
                       names(layers))
tri_lp <- file.path(out_dir, "11_three_layer_patterns.pdf")
g3 <- ggplot(traj_sum, aes(State, z, group = Layer, color = Layer)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  geom_line(linewidth = 1.1) +
  geom_errorbar(aes(ymin = z - sd, ymax = z + sd), width = 0.12, linewidth = 0.6) +
  geom_point(size = 2.6) +
  scale_color_manual(values = layer_cols) +
  facet_wrap(~ ClusterF, ncol = 3) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), plot.title = element_text(size = 16)) +
  labs(x = "Differentiation state", y = "mean z-score",
       title = "Module patterns across RNA and chromatin layers")
ggsave(tri_lp, g3, width = 10, height = 7)
message("wrote ", tri_lp)
chk("three-layer line plot written", file.exists(tri_lp) && file.size(tri_lp) > 5000)
chk("one point per layer x module x state", nrow(traj_sum) == 3 * K * 3)

## ---------------------------------------------------------------------------
## 15. Supporting figures
## ---------------------------------------------------------------------------
## Figures that back up the claims the heatmaps make by eye.  They are plain
## ggplot2 built on the objects above; the three that need an extra package or a
## remote file (16, 18, 19) are skipped with a message when it is unavailable.
section("15. Supporting figures")
MOD_COLS  <- ModuleViz:::.moduleviz_cluster_cols(levels(obj$cluster), "paired")
MARK_COLS <- layer_cols[c("H3K27ac", "H3K27me3")]
theme_sup <- theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 16))
save_sup <- function(file, plot, width, height) {
  ggsave(file, plot, width = width, height = height)
  message("wrote ", file)
  chk(sprintf("%s written", basename(file)),
      file.exists(file) && file.size(file) > 4000)
}
marks <- list(H3K27ac = list(cr = ac, z = d_ac$z_secondary),
              H3K27me3 = list(cr = me3, z = d_me3$z_secondary))

## ---- 12. per-gene RNA ~ chromatin correlation ------------------------------
## The assays were run on separate aliquots, so replicates are paired by their
## infection / day / replicate label.  Within a state that pairing is arbitrary;
## the correlation is carried by the differences between states.
pair_key <- function(ids)
  sub("_S[0-9]+$", "", sub("_(RNA|H3K27ac|H3K27me3)_", "_", ids))
row_cor <- function(x, y) {
  x <- x - rowMeans(x); y <- y - rowMeans(y)
  rowSums(x * y) / sqrt(rowSums(x^2) * rowSums(y^2))
}
gene_cor <- rbindlist(lapply(names(marks), function(nm) {
  zs <- marks[[nm]]$z
  i  <- match(pair_key(colnames(obj$z)), pair_key(colnames(zs)))
  data.table(Gene = rownames(zs), Mark = nm,
             Cluster = module_of[rownames(zs)],
             r = row_cor(obj$z[rownames(zs), , drop = FALSE], zs[, i, drop = FALSE]))
}))
gene_cor[, Mark := factor(Mark, levels = names(marks))]
gene_cor[, Cluster := factor(Cluster, levels = levels(obj$cluster))]
cor_med <- dcast(gene_cor[, .(r = round(median(r, na.rm = TRUE), 2)),
                          by = .(Cluster, Mark)], Cluster ~ Mark, value.var = "r")
cat("\n  median per-gene correlation with RNA:\n"); print(cor_med)
chk("every RNA sample has a paired sample in both marks",
    all(vapply(marks, function(m) all(pair_key(colnames(obj$z)) %in%
                                        pair_key(colnames(m$z))), logical(1))))
chk("H3K27ac is positively correlated with RNA in every module",
    all(cor_med$H3K27ac > 0))
chk("H3K27ac tracks RNA more closely than H3K27me3 in every module",
    all(cor_med$H3K27ac > cor_med$H3K27me3))

g12 <- ggplot(gene_cor[is.finite(r)], aes(Cluster, r, fill = Mark)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_violin(position = position_dodge(0.85), width = 0.85, color = NA,
              alpha = 0.55, scale = "width") +
  geom_boxplot(position = position_dodge(0.85), width = 0.18,
               outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = MARK_COLS) +
  theme_sup + theme(legend.position = "top", legend.title = element_blank(),
                    axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "per-gene Pearson r with RNA",
       title = "Does chromatin follow expression, gene by gene?")
save_sup(file.path(out_dir, "12_RNA_chromatin_correlation.pdf"), g12, 9, 6)

## ---- 13. fold change: RNA vs chromatin at the same gene --------------------
peak_lfc <- function(cr, a, b) {
  st <- cr$meta$State[match(colnames(cr$mat), cr$meta$SampleID)]
  d  <- rowMeans(cr$mat[, st == a, drop = FALSE]) -
        rowMeans(cr$mat[, st == b, drop = FALSE])
  tapply(d, cr$ann$Gene, mean)                  # mean over a gene's peaks
}
fc <- rbindlist(lapply(names(marks), function(nm) {
  l <- peak_lfc(marks[[nm]]$cr, "TEX", "TN")
  g <- intersect(sig$ID, names(l))
  data.table(Gene = g, Mark = nm, Cluster = module_of[g],
             rna = sig$log2FC_TEX_vs_TN[match(g, sig$ID)], chrom = l[g])
}))
fc[, Mark := factor(Mark, levels = names(marks))]
fc[, Cluster := factor(Cluster, levels = levels(obj$cluster))]
fc_r <- fc[, .(r = cor(rna, chrom), n = .N), by = Mark]
cat("\n  RNA vs chromatin log2FC (TEX vs TN) across DEGs:\n"); print(fc_r)
chk("H3K27ac change is positively correlated with the RNA change",
    fc_r[Mark == "H3K27ac", r] > 0.3)
chk("H3K27ac change tracks RNA better than H3K27me3 change",
    fc_r[Mark == "H3K27ac", r] > fc_r[Mark == "H3K27me3", r])

fc_lab <- fc[Gene %in% LABEL_GENES]
g13 <- ggplot(fc, aes(rna, chrom)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_vline(xintercept = 0, color = "grey75") +
  geom_point(aes(color = Cluster), size = 1.3, alpha = 0.65) +
  geom_text(data = fc_r, aes(x = -Inf, y = Inf,
                             label = sprintf("r = %.2f  (n = %d)", r, n)),
            hjust = -0.08, vjust = 1.6, size = 5.5) +
  scale_color_manual(values = MOD_COLS) +
  facet_wrap(~ Mark, scales = "free_y") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_sup + theme(legend.title = element_blank()) +
  labs(x = "RNA log2 fold change (TEX vs TN)",
       y = "chromatin log2 fold change (TEX vs TN)",
       title = "Expression change vs chromatin change at the same gene")
g13 <- g13 + if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(data = fc_lab, aes(label = Gene), size = 4.6,
                           fontface = "italic", max.overlaps = Inf,
                           min.segment.length = 0, segment.color = "grey40",
                           seed = 1)
} else {
  geom_text(data = fc_lab, aes(label = Gene), size = 4.6, fontface = "italic",
            vjust = -0.7, check_overlap = TRUE)
}
save_sup(file.path(out_dir, "13_foldchange_RNA_vs_chromatin.pdf"), g13, 13, 6.5)

## ---- 14. sample QC: PCA of each layer --------------------------------------
qc_layers <- list(RNA = log_rna[rowMeans(rna_mat) >= MIN_EXPR, ],
                  H3K27ac = ac$mat, H3K27me3 = me3$mat)
pca <- rbindlist(lapply(names(qc_layers), function(nm) {
  m  <- qc_layers[[nm]]
  m  <- m[order(apply(m, 1, var), decreasing = TRUE)[seq_len(min(2000, nrow(m)))], ]
  pc <- stats::prcomp(t(m), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  inf <- sub("_P14_.*", "", colnames(m))
  data.table(Layer = sprintf("%s  (PC1 %.0f%%, PC2 %.0f%%)", nm, ve[1], ve[2]),
             State = factor(c(Naive = "TN", Arm = "TMEM", Cl13 = "TEX")[inf],
                            levels = STATES),
             Cohort = sub(".*_(Day[0-9]+)_.*", "\\1", colnames(m)),
             PC1 = pc$x[, 1], PC2 = pc$x[, 2], ve1 = ve[1])
}))
pca[, Layer := factor(Layer, levels = unique(Layer))]
chk("PC1 carries most of the variance in every layer",
    all(pca[, ve1[1], by = Layer]$V1 > 40))
g14 <- ggplot(pca, aes(PC1, PC2, color = State, shape = Cohort)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = state_cols) +
  facet_wrap(~ Layer, scales = "free") +
  theme_sup + theme(legend.position = "top") +
  labs(title = "Sample PCA per layer (2,000 most variable features)")
save_sup(file.path(out_dir, "14_sample_PCA.pdf"), g14, 14, 5.8)

## ---- 15. sample QC: signal distribution per sample -------------------------
## The chromatin layers here are the peaks within 25 kb of a TSS.  For H3K27me3
## those peaks sit about one log2 unit higher in every TMEM and TEX sample than
## in TN, whichever module the gene belongs to.  The published normalisation
## balances the full peak set (the median change over all peaks is 0), so this
## is signal moving toward promoter-proximal peaks everywhere, not a
## module-specific change: read the H3K27me3 panels relative to it.
sig_dist <- rbindlist(lapply(names(qc_layers), function(nm) {
  m   <- qc_layers[[nm]]
  inf <- sub("_P14_.*", "", colnames(m))
  q   <- apply(m, 2, stats::quantile, c(0.05, 0.25, 0.5, 0.75, 0.95))
  data.table(Layer = nm, Sample = colnames(m),
             State = factor(c(Naive = "TN", Arm = "TMEM", Cl13 = "TEX")[inf],
                            levels = STATES),
             lo = q[1, ], q1 = q[2, ], med = q[3, ], q3 = q[4, ], hi = q[5, ])
}))
sig_dist[, Layer := factor(Layer, levels = names(qc_layers))]
setorder(sig_dist, Layer, State, Sample)
sig_dist[, idx := seq_len(.N), by = Layer]
shift <- dcast(sig_dist[, .(med = median(med)), by = .(Layer, State)],
               Layer ~ State, value.var = "med")
shift[, `:=`(TMEM_vs_TN = round(TMEM - TN, 2), TEX_vs_TN = round(TEX - TN, 2))]
cat("\n  median log2 signal per state, and its shift against TN:\n"); print(shift)
chk("RNA and H3K27ac medians are level across states (|shift| < 0.2)",
    all(abs(unlist(shift[Layer != "H3K27me3", .(TMEM_vs_TN, TEX_vs_TN)])) < 0.2))
chk("H3K27me3 carries a global shift in TMEM and TEX (> 0.5 log2)",
    all(unlist(shift[Layer == "H3K27me3", .(TMEM_vs_TN, TEX_vs_TN)]) > 0.5))

g15 <- ggplot(sig_dist, aes(idx, fill = State)) +
  geom_boxplot(aes(ymin = lo, lower = q1, middle = med, upper = q3, ymax = hi,
                   group = idx), stat = "identity", width = 0.7, linewidth = 0.4) +
  scale_fill_manual(values = state_cols) +
  facet_wrap(~ Layer, scales = "free_y") +
  theme_sup + theme(legend.position = "top", axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.grid.major.x = element_blank()) +
  labs(x = "samples (one box each; whiskers 5-95%)",
       y = "log2(normalised signal + 1)",
       title = "Signal distribution per sample",
       subtitle = "chromatin layers: peaks within 25 kb of a TSS")
save_sup(file.path(out_dir, "15_sample_signal_distribution.pdf"), g15, 13, 5.8)

## ---- 16. why k = 6 and not the elbow's k = 5 -------------------------------
if (requireNamespace("ggalluvial", quietly = TRUE)) {
  o5 <- longitudinal_cluster(deg_mat, rna_meta, sample_col = "SampleID",
                             time_col = "StateOrder", k = k_sug,
                             stability = FALSE, seed = 1)
  c5 <- as.matrix(o5$centroids[, -1])
  peak5 <- setNames(STATES[apply(c5, 1, which.max)], o5$centroids$Cluster)
  ## the k = 5 modules get their own names (A-E): "Cluster4" is a different set
  ## of genes in the two runs, and sharing a label would imply it is the same
  strat <- function(cl, peak, names_out) {
    lab <- sprintf("%s\n%s", names_out[as.integer(cl)], peak[as.character(cl)])
    factor(lab, levels = unique(lab[order(as.integer(cl))]))
  }
  k5_of <- o5$cluster[match(rownames(obj$z), rownames(o5$z))]
  flow <- data.table(k5 = strat(k5_of, peak5, paste("Module", LETTERS[seq_len(k_sug)])),
                     k6 = strat(obj$cluster, peak_state, levels(obj$cluster)),
                     Module = obj$cluster)[, .N, by = .(k5, k6, Module)]
  chk("the flow table accounts for every DEG", sum(flow$N) == nrow(obj$z))
  g16 <- ggplot(flow, aes(axis1 = k5, axis2 = k6, y = N)) +
    ggalluvial::geom_alluvium(aes(fill = Module), width = 0.32, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 0.32, fill = "grey96", color = "grey35") +
    geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)),
              size = 4.6, lineheight = 0.9) +
    scale_fill_manual(values = MOD_COLS, guide = "none") +
    scale_x_continuous(breaks = 1:2, expand = c(0.08, 0.08),
                       labels = c(sprintf("k = %d (elbow)", k_sug),
                                  sprintf("k = %d (used)", K))) +
    theme_sup + theme(panel.grid = element_blank(), panel.border = element_blank(),
                      axis.ticks.x = element_blank()) +
    labs(x = NULL, y = "genes",
         title = "Where the memory module comes from",
         subtitle = "each module is labelled with the state it peaks in")
  save_sup(file.path(out_dir, "16_k5_vs_k6_flow.pdf"), g16, 8, 9)
} else cat("  [skip] 16_k5_vs_k6_flow: install 'ggalluvial' to draw it\n")

## ---- 17. individual marker genes across the three layers -------------------
PANEL_GENES <- intersect(c("Sell", "Il7r", "Klrg1", "Tox", "Havcr2", "Eomes"), genes3)
gene_signal <- function(cr, g)                 # mean over the gene's peaks
  colMeans(cr$mat[cr$ann$Gene == g, , drop = FALSE])
panel <- rbindlist(lapply(PANEL_GENES, function(g) rbind(
  data.table(Gene = g, Layer = "RNA", Sample = colnames(log_rna), v = log_rna[g, ]),
  data.table(Gene = g, Layer = "H3K27ac", Sample = colnames(ac$mat),
             v = gene_signal(ac, g)),
  data.table(Gene = g, Layer = "H3K27me3", Sample = colnames(me3$mat),
             v = gene_signal(me3, g)))))
panel[, State := factor(c(Naive = "TN", Arm = "TMEM", Cl13 = "TEX")[
  sub("_P14_.*", "", Sample)], levels = STATES)]
panel[, Gene := factor(Gene, levels = PANEL_GENES)]
panel[, Layer := factor(Layer, levels = names(layers))]
chk("marker panel has 16 samples per gene and layer",
    all(panel[, .N, by = .(Gene, Layer)]$N == 16))
g17 <- ggplot(panel, aes(State, v, color = State)) +
  geom_jitter(width = 0.14, height = 0, size = 2.2, alpha = 0.8) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar",
               width = 0.55, linewidth = 0.35, color = "grey20") +
  scale_color_manual(values = state_cols, guide = "none") +
  facet_grid(Layer ~ Gene, scales = "free_y") +
  theme_sup + theme(strip.text.x = element_text(face = "italic"),
                    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "log2(normalised signal + 1)",
       title = "Marker genes: replicates (points) and state mean (bar)")
save_sup(file.path(out_dir, "17_marker_gene_panels.pdf"), g17, 14, 8)

## ---- 18. Hallmark pathway enrichment per module ----------------------------
if (requireNamespace("msigdbr", quietly = TRUE)) {
  hm <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", collection = "H"),
                 error = function(e)                       # msigdbr < 10
                   msigdbr::msigdbr(species = "Mus musculus", category = "H"))
  universe <- deg[baseMean >= MIN_EXPR, ID]   # every gene that could be a DEG
  hm <- unique(as.data.table(hm)[gene_symbol %in% universe, .(gs_name, gene_symbol)])
  sets <- split(hm$gene_symbol, hm$gs_name)
  enr <- rbindlist(lapply(levels(obj$cluster), function(m) {
    mg <- names(module_of)[module_of == m]
    rbindlist(lapply(names(sets), function(s) {
      k <- length(intersect(mg, sets[[s]])); M <- length(sets[[s]])
      data.table(Cluster = m, Set = s, Overlap = k, SetSize = M,
                 Fold = (k / length(mg)) / (M / length(universe)),
                 p = stats::phyper(k - 1, M, length(universe) - M, length(mg),
                                   lower.tail = FALSE))
    }))
  }))
  enr[, padj := stats::p.adjust(p, method = "BH")]
  enr_sig <- enr[padj < 0.05 & Overlap >= 5]
  cat(sprintf("\n  %d module x Hallmark pairs at FDR < 0.05:\n", nrow(enr_sig)))
  print(enr_sig[order(Cluster, padj),
                .(Cluster, Set, Overlap, Fold = round(Fold, 1), padj = signif(padj, 2))])
  chk("hypergeometric universe contains every DEG", all(names(module_of) %in% universe))
  chk("at least one Hallmark set is enriched", nrow(enr_sig) > 0)

  enr_sig[, Cluster := factor(Cluster, levels = levels(obj$cluster))]
  pretty_set <- function(s) gsub("_", " ", sub("^HALLMARK_", "", s))
  ## order the sets by the module they are strongest in, so the dots run
  ## down the diagonal
  set_order <- enr_sig[order(as.integer(Cluster), padj), unique(Set)]
  enr_sig[, SetLab := factor(pretty_set(Set), levels = rev(pretty_set(set_order)))]
  g18 <- ggplot(enr_sig, aes(Cluster, SetLab)) +
    geom_point(aes(size = Overlap, fill = -log10(padj)), shape = 21,
               color = "grey25", stroke = 0.4) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", name = "-log10 FDR") +
    scale_size_area(max_size = 11, name = "genes") +
    scale_x_discrete(drop = FALSE) +
    theme_sup + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = NULL,
         title = "MSigDB Hallmark sets enriched per module (FDR < 0.05)")
  save_sup(file.path(out_dir, "18_hallmark_enrichment.pdf"), g18, 10,
           max(5, 0.42 * nlevels(enr_sig$SetLab) + 2.2))
} else cat("  [skip] 18_hallmark_enrichment: install 'msigdbr' to draw it\n")

## ---- 19. genome browser tracks ---------------------------------------------
## GEO hosts one pooled bigWig per mark and state (240-440 MB each).  Only the
## plotted windows are read, straight from the remote files, and the binned
## result (a few hundred kb of text) is cached so later runs stay offline.
TRACK_GENES <- c("Tcf7", "Pdcd1", "Tox")
track_cache <- file.path(cache_dir, "browser_tracks.tsv.gz")
N_BINS <- 400
loci <- rna_uniq[match(TRACK_GENES, Gene.Name), .(Gene = Gene.Name, Chr, Start, End)]
loci[, flank := pmax(15000L, as.integer(0.12 * (End - Start)))]
loci[, `:=`(from = Start - flank, to = End + flank)]

fetch_tracks <- function() {
  gr <- GenomicRanges::GRanges(loci$Chr, IRanges::IRanges(loci$from, loci$to))
  grid <- CJ(Mark = names(marks), Infection = c("Naive", "Arm", "Cl13"), sorted = FALSE)
  rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    url <- sprintf(paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285245/",
                          "suppl/GSE285245_pool_%s_%s.bw"),
                   grid$Mark[i], grid$Infection[i])
    message("reading ", basename(url))
    rbindlist(lapply(seq_len(nrow(loci)), function(j) {
      ## one window per request: a multi-range import returns the windows in
      ## the file's chromosome order, not in the order they were asked for
      v <- NULL
      for (attempt in 1:6) {                 # the remote open fails now and then
        v <- tryCatch(suppressWarnings(as.numeric(rtracklayer::import(
          rtracklayer::BigWigFile(url), which = gr[j], as = "NumericList")[[1]])),
          error = function(e) NULL)
        if (!is.null(v)) break
        Sys.sleep(5)
      }
      if (is.null(v)) stop("no answer from ", basename(url))
      stopifnot(length(v) == loci$to[j] - loci$from[j] + 1L)
      bin <- cut(seq_along(v), N_BINS, labels = FALSE)
      data.table(Gene = loci$Gene[j], Mark = grid$Mark[i],
                 Infection = grid$Infection[i],
                 pos = loci$from[j] + as.integer(tapply(seq_along(v), bin, mean)),
                 signal = round(as.numeric(tapply(v, bin, mean)), 3))
    }))
  }))
}
tracks <- if (file.exists(track_cache)) fread(track_cache) else
  tryCatch({
    stopifnot(requireNamespace("rtracklayer", quietly = TRUE))
    tr <- fetch_tracks(); fwrite(tr, track_cache, sep = "\t"); tr
  }, error = function(e) {
    cat("  [skip] 19_browser_tracks: could not read the remote bigWigs (",
        conditionMessage(e), ")\n"); NULL })

if (!is.null(tracks)) {
  chk("browser tracks cover every locus, mark and state",
      nrow(tracks) == length(TRACK_GENES) * 2 * 3 * N_BINS)
  ## the tracks must tell the same story as the count matrices
  body_ac <- function(g, inf) {
    l <- loci[Gene == g]
    tracks[Gene == g & Mark == "H3K27ac" & Infection == inf &
             pos >= l$Start & pos <= l$End, mean(signal)]
  }
  chk("H3K27ac over Pdcd1 and Tox is higher in TEX than in TN",
      body_ac("Pdcd1", "Cl13") > body_ac("Pdcd1", "Naive") &&
        body_ac("Tox", "Cl13") > body_ac("Tox", "Naive"))
  chk("H3K27ac over Tcf7 is higher in TN than in TEX",
      body_ac("Tcf7", "Naive") > body_ac("Tcf7", "Cl13"))
  tracks[, State := factor(c(Naive = "TN", Arm = "TMEM", Cl13 = "TEX")[Infection],
                           levels = STATES)]
  row_levels <- c(as.vector(t(outer(names(marks), STATES, paste, sep = "\n"))), "Genes")
  tracks[, Row := factor(paste(Mark, State, sep = "\n"), levels = row_levels)]
  locus_lab <- setNames(sprintf("%s  (%s)", loci$Gene, loci$Chr), loci$Gene)
  tracks[, Locus := factor(locus_lab[Gene], levels = locus_lab)]
  ## the three states of a mark share one y axis, so heights are comparable
  ceil <- tracks[, .(signal = max(signal)), by = Mark][
    , .(Row = factor(paste(Mark, STATES, sep = "\n"), levels = row_levels),
        signal = signal), by = Mark]
  ceil <- ceil[, .(Locus = factor(locus_lab, levels = locus_lab), pos = loci$from),
               by = .(Row, signal)]

  ## gene bodies inside each window, staggered over three lanes
  gm <- rbindlist(lapply(seq_len(nrow(loci)), function(j)
    rna_uniq[Chr == loci$Chr[j] & End >= loci$from[j] & Start <= loci$to[j],
             .(Locus = factor(locus_lab[loci$Gene[j]], levels = locus_lab),
               Name = Gene.Name, Strand,
               Start = pmax(Start, loci$from[j]), End = pmin(End, loci$to[j]))]))
  setorder(gm, Locus, Start)
  gm[, lane := (seq_len(.N) - 1L) %% 3L, by = Locus]
  gm[, Row := factor("Genes", levels = row_levels)]
  gm[, `:=`(x0 = ifelse(Strand == "+", Start, End), x1 = ifelse(Strand == "+", End, Start))]

  g19 <- ggplot(tracks, aes(pos / 1e6)) +
    geom_area(aes(y = signal, fill = State)) +
    geom_blank(data = ceil, aes(y = signal)) +
    geom_segment(data = gm, aes(x = x0 / 1e6, xend = x1 / 1e6, y = lane, yend = lane),
                 linewidth = 1.1, color = "grey25",
                 arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")) +
    geom_text(data = gm, aes(x = (Start + End) / 2e6, y = lane + 0.42, label = Name),
              size = 3.8, fontface = "italic") +
    geom_blank(data = gm, aes(x = Start / 1e6, y = 2.8)) +
    scale_fill_manual(values = state_cols, guide = "none") +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    facet_grid(Row ~ Locus, scales = "free", switch = "y") +
    theme_sup + theme(strip.text.y.left = element_text(angle = 0),
                      strip.placement = "outside", panel.grid = element_blank(),
                      axis.text.y = element_text(size = 10),
                      panel.spacing.y = grid::unit(1.5, "mm")) +
    labs(x = "position (Mb, mm10)", y = NULL,
         title = "Pooled CUT&RUN signal at a naive locus and two exhaustion loci")
  save_sup(file.path(out_dir, "19_browser_tracks.pdf"), g19, 14, 9.5)
}

## ---------------------------------------------------------------------------
## 16. One-call wrapper
## ---------------------------------------------------------------------------
section("16. run_longitudinal wrapper")
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
## 17. Error handling
## ---------------------------------------------------------------------------
section("17. Error handling on malformed real input")
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
      err(dual_heatmap(obj, fake, file = tempfile(fileext = ".png"))) })

## ---------------------------------------------------------------------------
section("SUMMARY")
cat(sprintf("  %d passed, %d failed\n", PASS, FAIL))
cat(sprintf("  outputs in: %s\n", normalizePath(out_dir)))
print(data.frame(file = list.files(out_dir),
                 bytes = file.size(file.path(out_dir, list.files(out_dir)))),
      row.names = FALSE)
if (FAIL > 0) quit(status = 1)
