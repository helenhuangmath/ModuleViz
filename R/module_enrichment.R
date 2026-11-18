## =============================================================================
## module_enrichment.R
##
## Over-representation analysis of gene sets within temporal modules, plus the
## helpers that print the top terms beside each module of a heatmap.
##
## The test is a plain hypergeometric test against user-supplied gene sets, so
## it needs no annotation package and works for any organism or any custom
## collection (pathways, TF targets, marker lists, ...).
## =============================================================================

#' Normalise gene sets to a named list of character vectors.
#'
#' Accepts a named list, a two-column data.frame (term, gene), or the path to a
#' GMT file.
#' @noRd
.moduleviz_gene_sets <- function(gene_sets) {
    if (is.character(gene_sets) && length(gene_sets) == 1 &&
        is.null(names(gene_sets))) {
        if (!file.exists(gene_sets)) {
            stop("Gene-set file not found: ", gene_sets, call. = FALSE)
        }
        lines <- strsplit(readLines(gene_sets, warn = FALSE), "\t")
        lines <- lines[lengths(lines) >= 3]
        sets <- lapply(lines, function(x) x[-c(1, 2)])
        names(sets) <- vapply(lines, `[`, character(1), 1)
        gene_sets <- sets
    } else if (is.data.frame(gene_sets)) {
        if (ncol(gene_sets) < 2) {
            stop("A gene-set data frame needs two columns: term and gene.",
                 call. = FALSE)
        }
        gene_sets <- split(as.character(gene_sets[[2]]),
                           as.character(gene_sets[[1]]))
    }
    if (!is.list(gene_sets) || is.null(names(gene_sets)) ||
        any(!nzchar(names(gene_sets)))) {
        stop("`gene_sets` must be a named list of gene vectors, a two-column ",
             "data frame (term, gene), or the path to a GMT file.",
             call. = FALSE)
    }
    lapply(gene_sets, function(g) {
        g <- unique(as.character(g))
        g[!is.na(g) & nzchar(g)]
    })
}

#' Gene-set enrichment of every temporal module.
#'
#' Tests each module for over-representation of each gene set with a
#' one-sided hypergeometric test.  The result can be passed to the
#' `enrichment` argument of [longitudinal_heatmap()] and [multi_heatmap()] to
#' print the top terms beside each module, or drawn on its own with
#' [enrichment_dotplot()].
#'
#' @param obj a `longi` or `multi_longi` object.
#' @param gene_sets gene sets to test: a named list of gene vectors, a
#'   two-column data frame (term, gene), or the path to a GMT file.
#' @param universe background genes.  Defaults to every gene in `obj`, which
#'   is appropriate when the object holds all tested genes; supply the full
#'   set of measured genes when `obj` was built from a pre-selected subset.
#' @param gene_of optional function mapping feature IDs to gene symbols, for
#'   objects whose rows are not genes (e.g. `function(id) sub("^.*_", "", id)`
#'   for `peak12_Cd8a`).  Ignored for `multi_longi` objects, which already
#'   carry gene symbols.
#' @param min_size,max_size gene sets with fewer or more genes in the universe
#'   are skipped.
#' @param min_overlap module / gene-set pairs sharing fewer genes are not
#'   tested.
#' @param p_adjust_method multiple-testing correction passed to
#'   [stats::p.adjust()], applied across all tested module / gene-set pairs.
#' @return A `data.table` with one row per tested module / gene-set pair:
#'   `Cluster`, `Term`, `Overlap`, `ModuleSize`, `SetSize`, `Universe`,
#'   `FoldEnrichment`, `PValue`, `Padj`, and `Genes` (the shared genes),
#'   ordered by module and p-value.
#' @examples
#' data("heatmap_example", package = "ModuleViz")
#' obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
#'                             time_col = "TimeHr", group_col = "Condition",
#'                             k = 4, stability = FALSE)
#' genes <- rownames(heatmap_example$mat)
#' sets <- list(Early_response = genes[1:15], Late_response = genes[31:45],
#'              Housekeeping = genes[seq(2, 60, by = 4)])
#' enr <- module_enrichment(obj, sets, min_size = 3)
#' head(enr)
module_enrichment <- function(obj, gene_sets, universe = NULL, gene_of = NULL,
                              min_size = 5, max_size = 500, min_overlap = 2,
                              p_adjust_method = "BH") {
    if (!inherits(obj, c("longi", "multi_longi"))) {
        stop("`obj` must be a 'longi' or 'multi_longi' object.", call. = FALSE)
    }
    gene_sets <- .moduleviz_gene_sets(gene_sets)
    if (inherits(obj, "multi_longi")) {
        genes <- as.character(obj$genes)
    } else {
        genes <- rownames(obj$z)
        if (!is.null(gene_of)) {
            genes <- as.character(gene_of(genes))
        }
    }
    cluster <- obj$cluster
    ok <- !is.na(genes) & nzchar(genes)
    modules <- lapply(split(genes[ok], cluster[ok]), unique)

    if (is.null(universe)) {
        universe <- unique(genes[ok])
    } else {
        universe <- unique(as.character(universe))
        outside <- setdiff(unique(genes[ok]), universe)
        if (length(outside)) {
            message(length(outside), " module gene(s) are not in `universe` ",
                    "and are ignored.")
        }
        modules <- lapply(modules, intersect, universe)
    }
    n_universe <- length(universe)

    gene_sets <- lapply(gene_sets, intersect, universe)
    set_size <- lengths(gene_sets)
    gene_sets <- gene_sets[set_size >= min_size & set_size <= max_size]
    empty <- data.table(Cluster = character(0), Term = character(0),
                        Overlap = integer(0), ModuleSize = integer(0),
                        SetSize = integer(0), Universe = integer(0),
                        FoldEnrichment = numeric(0), PValue = numeric(0),
                        Padj = numeric(0), Genes = character(0))
    if (!length(gene_sets)) {
        warning("No gene set has between ", min_size, " and ", max_size,
                " genes in the universe; nothing was tested.", call. = FALSE)
        return(empty)
    }

    res <- rbindlist(lapply(names(modules), function(cl) {
        mod <- modules[[cl]]
        rbindlist(lapply(names(gene_sets), function(term) {
            set <- gene_sets[[term]]
            shared <- intersect(mod, set)
            if (length(shared) < min_overlap) {
                return(NULL)
            }
            p <- stats::phyper(length(shared) - 1, length(set),
                               n_universe - length(set), length(mod),
                               lower.tail = FALSE)
            expected <- length(mod) * length(set) / n_universe
            data.table(Cluster = cl, Term = term, Overlap = length(shared),
                       ModuleSize = length(mod), SetSize = length(set),
                       Universe = n_universe,
                       FoldEnrichment = round(length(shared) / expected, 3),
                       PValue = p, Genes = paste(shared, collapse = ","))
        }))
    }))
    if (!nrow(res)) {
        return(empty)
    }
    res[, Padj := stats::p.adjust(PValue, method = p_adjust_method)]
    res[, Cluster := factor(Cluster, levels = levels(cluster))]
    setorder(res, Cluster, PValue)
    res[, Cluster := as.character(Cluster)]
    setcolorder(res, names(empty))
    res[]
}

#' Check an enrichment table and return it as a data.table.
#' @noRd
.moduleviz_check_enrichment <- function(enrichment) {
    enrichment <- as.data.table(enrichment)
    need <- c("Cluster", "Term")
    absent <- setdiff(need, names(enrichment))
    if (length(absent)) {
        stop("`enrichment` must have the columns Cluster and Term (and ",
             "ideally Padj); missing: ", paste(absent, collapse = ", "),
             call. = FALSE)
    }
    if (!"Padj" %in% names(enrichment)) {
        enrichment[, Padj := if ("PValue" %in% names(enrichment)) PValue
                             else NA_real_]
    }
    enrichment
}

#' Top enriched terms per module, as one text block per module.
#'
#' @return character vector named by module; "" where nothing passes.
#' @noRd
.moduleviz_enrichment_labels <- function(enrichment, modules, top_n = 3,
                                         padj_cutoff = 0.05, max_chars = 45) {
    enrichment <- .moduleviz_check_enrichment(enrichment)
    keep <- enrichment[is.na(Padj) | Padj <= padj_cutoff]
    setorder(keep, Padj, na.last = TRUE)
    vapply(modules, function(cl) {
        terms <- head(keep[Cluster == cl]$Term, top_n)
        if (!length(terms)) {
            return("")
        }
        long <- nchar(terms) > max_chars
        terms[long] <- paste0(substr(terms[long], 1, max_chars - 3), "...")
        paste(terms, collapse = "\n")
    }, character(1))
}

#' Right-hand heatmap annotation listing the top enriched terms per module.
#'
#' @return a ComplexHeatmap annotation function, or NULL when no module has a
#'   term to show.
#' @noRd
.moduleviz_enrichment_anno <- function(enrichment, modules, top_n = 3,
                                       padj_cutoff = 0.05, fontsize = 9) {
    labels <- .moduleviz_enrichment_labels(enrichment, modules, top_n,
                                           padj_cutoff)
    if (!any(nzchar(labels))) {
        message("No enriched term passes `enrichment_padj`; ",
                "the enrichment annotation is omitted.")
        return(NULL)
    }
    gp <- gpar(fontsize = fontsize)
    text_lines <- unlist(strsplit(labels[nzchar(labels)], "\n"))
    anno_block(which = "row", gp = gpar(fill = NA, col = NA),
               labels = unname(labels),
               labels_gp = gp, labels_rot = 0, labels_just = "left",
               labels_offset = unit(1.5, "mm"),
               width = max_text_width(text_lines, gp = gp) + unit(4, "mm"))
}

#' Dot plot of the top enriched gene sets per module.
#'
#' @param enrichment the table returned by [module_enrichment()] (or any data
#'   frame with `Cluster`, `Term`, `Padj`, and `FoldEnrichment` columns).
#' @param top_n number of terms shown per module, ranked by adjusted p-value.
#' @param padj_cutoff terms with a larger adjusted p-value are not shown.
#' @param base_size ggplot base font size.
#' @param title optional plot title; none by default.
#' @param file optional output path (any [ggplot2::ggsave()] format).
#' @param width,height size in inches; auto-scaled when `NULL`.
#' @return a `ggplot` object (also saved if `file` is given).
#' @examples
#' data("heatmap_example", package = "ModuleViz")
#' obj <- longitudinal_cluster(heatmap_example$mat, heatmap_example$meta,
#'                             time_col = "TimeHr", group_col = "Condition",
#'                             k = 4, stability = FALSE)
#' genes <- rownames(heatmap_example$mat)
#' sets <- list(Early_response = genes[1:15], Late_response = genes[31:45])
#' enr <- module_enrichment(obj, sets, min_size = 3)
#' enrichment_dotplot(enr, padj_cutoff = 1)
enrichment_dotplot <- function(enrichment, top_n = 5, padj_cutoff = 0.05,
                               base_size = 10,
                               title = NULL,
                               file = NULL, width = NULL, height = NULL) {
    enrichment <- .moduleviz_check_enrichment(enrichment)
    if (!"FoldEnrichment" %in% names(enrichment)) {
        enrichment[, FoldEnrichment := 1]
    }
    keep <- enrichment[!is.na(Padj) & Padj <= padj_cutoff]
    if (!nrow(keep)) {
        stop("No term passes `padj_cutoff` = ", padj_cutoff, ".",
             call. = FALSE)
    }
    setorder(keep, Padj)
    keep <- keep[, head(.SD, top_n), by = Cluster]
    modules <- unique(enrichment$Cluster)
    keep[, Cluster := factor(Cluster, levels = modules)]
    keep[, Term := factor(Term, levels = rev(unique(Term[order(Cluster)])))]
    keep[, Significance := -log10(pmax(Padj, .Machine$double.xmin))]

    gg <- ggplot(keep, aes(Cluster, Term)) +
        geom_point(aes(size = FoldEnrichment, color = Significance)) +
        scale_color_gradient(low = "#9ecae1", high = "#08306b") +
        scale_size_continuous(range = c(1.5, 4.5)) +
        scale_x_discrete(drop = FALSE) +
        .moduleviz_theme(base_size) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.title = element_text(size = base_size * 0.9),
              legend.key.height = unit(3, "mm"),
              legend.box = "horizontal") +
        labs(x = NULL, y = NULL, size = "Fold",
             color = "-log10 padj", title = title)

    if (!is.null(file)) {
        if (is.null(width)) {
            width <- 2.4 + 0.4 * length(modules) +
                0.07 * max(nchar(as.character(keep$Term)))
        }
        if (is.null(height)) {
            height <- max(2.2, 1.2 + 0.22 * nlevels(keep$Term))
        }
        ggsave(file, gg, width = width, height = height, limitsize = FALSE)
        message("wrote ", file)
    }
    gg
}
