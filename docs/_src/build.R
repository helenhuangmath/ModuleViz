#!/usr/bin/env Rscript
## Build the ModuleViz documentation site into docs/.
##
## Run from the package root:
##     Rscript docs/_src/build.R
##
## Sources:
##   README.md            sliced by "## " heading into the guide pages
##   NEWS.md              release notes
##   man/*.Rd             function reference (tools::Rd2HTML)
##   docs/_src/pages/*.md hand-written pages (home, install, formats, Q & A,
##                        vignettes)
## Needs pandoc on the PATH. Output: docs/*.html, docs/reference/*.html,
## docs/search-index.js. docs/assets/ and docs/articles/ are left untouched.

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "DESCRIPTION")))
src <- file.path(root, "docs", "_src")
out <- file.path(root, "docs")
dir.create(file.path(out, "reference"), showWarnings = FALSE)

desc <- read.dcf(file.path(root, "DESCRIPTION"))
version <- desc[1, "Version"]
repo <- "https://github.com/helenhuangmath/ModuleViz"

## ---------------------------------------------------------------------------
## Site map. `readme` names the README "## " sections a page is built from;
## `file` names a hand-written page; `news` builds the release notes.
## ---------------------------------------------------------------------------
pages <- list(
    list(id = "index", title = "ModuleViz", file = "index.md"),
    list(id = "install", title = "Install ModuleViz", file = "install.md"),
    list(id = "tutorial", title = "Tutorial",
         readme = c("Quick Example: Cluster Temporal Modules",
                    "Longitudinal Module Analysis"),
         keep_h2 = TRUE),
    list(id = "functions", title = "Function guide",
         readme = "Function Tutorials"),
    list(id = "multiomics", title = "Multi-omics analysis",
         readme = "Multi-Omics Analysis"),
    list(id = "coherence", title = "Module coherence",
         readme = "Module Coherence"),
    list(id = "singlecell", title = "Single-cell data",
         readme = "Single-Cell Data"),
    list(id = "cd8", title = "Published data example",
         readme = "Published Data Example"),
    list(id = "public-data", title = "Public time-course data",
         readme = "Public Time-Course Data"),
    list(id = "formats", title = "Input and output formats",
         file = "formats.md"),
    list(id = "qa", title = "Common Q & A", file = "qa.md"),
    list(id = "vignettes", title = "Vignettes", file = "vignettes.md"),
    list(id = "reference/index", title = "Function reference",
         reference_index = TRUE),
    list(id = "changelog", title = "Release notes", news = TRUE)
)

## README anchors that now live on another page.
anchor_map <- c(
    "#multi-omics-analysis" = "multiomics.html",
    "#single-cell-data" = "singlecell.html",
    "#module-coherence" = "coherence.html",
    "#installation" = "install.html"
)

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
}

strip_tags <- function(x) {
    x <- gsub("(?s)<(script|style)[^>]*>.*?</\\1>", " ", x, perl = TRUE)
    x <- gsub("<[^>]+>", " ", x)
    x <- gsub("&lt;", "<", x, fixed = TRUE)
    x <- gsub("&gt;", ">", x, fixed = TRUE)
    x <- gsub("&quot;", "\"", x, fixed = TRUE)
    x <- gsub("&#39;", "'", x, fixed = TRUE)
    x <- gsub("&amp;", "&", x, fixed = TRUE)
    trimws(gsub("[[:space:]]+", " ", x))
}

js_string <- function(x) {
    x <- gsub("\\", "\\\\", x, fixed = TRUE)
    x <- gsub("\"", "\\\"", x, fixed = TRUE)
    x <- gsub("</", "<\\/", x, fixed = TRUE)
    x <- gsub("[\r\n\t]+", " ", x)
    paste0("\"", x, "\"")
}

## Lines of README.md / NEWS.md that sit inside a fenced code block.
in_fence <- function(lines) {
    cumsum(grepl("^```", lines)) %% 2 == 1 | grepl("^```", lines)
}

## Body of one "## <name>" README section (without its heading line).
readme_lines <- readLines(file.path(root, "README.md"), warn = FALSE)
readme_fenced <- in_fence(readme_lines)
readme_section <- function(name) {
    h2 <- which(grepl("^## ", readme_lines) & !readme_fenced)
    start <- h2[sub("^## +", "", readme_lines[h2]) == name]
    if (length(start) != 1L)
        stop("README section not found (or not unique): ", name)
    later <- h2[h2 > start]
    end <- if (length(later)) later[1] - 1L else length(readme_lines)
    readme_lines[seq(start + 1L, end)]
}

## Move headings up `by` levels outside code fences.
promote <- function(lines, by = 1L) {
    fenced <- in_fence(lines)
    hit <- grepl("^#{2,6} ", lines) & !fenced
    lines[hit] <- substring(lines[hit], by + 1L)
    lines
}

fix_paths <- function(lines) {
    lines <- gsub("](inst/extdata/results/", "](assets/results/", lines,
                  fixed = TRUE)
    lines <- gsub("](vignettes/figures/", "](assets/figures/", lines,
                  fixed = TRUE)
    for (a in names(anchor_map))
        lines <- gsub(paste0("](", a, ")"), paste0("](", anchor_map[[a]], ")"),
                      lines, fixed = TRUE)
    lines
}

pandoc_template <- file.path(src, "fragment.template")
pandoc <- function(md_lines) {
    infile <- tempfile(fileext = ".md")
    outfile <- tempfile(fileext = ".html")
    writeLines(md_lines, infile, useBytes = TRUE)
    status <- system2("pandoc", c(
        "--from", "gfm", "--to", "html5", "--toc", "--toc-depth=3",
        "--wrap=none", "--template", shQuote(pandoc_template),
        "--output", shQuote(outfile), shQuote(infile)))
    if (status != 0) stop("pandoc failed")
    txt <- paste(readLines(outfile, warn = FALSE, encoding = "UTF-8"),
                 collapse = "\n")
    parts <- strsplit(txt, "<!--MV-SPLIT-->", fixed = TRUE)[[1]]
    toc <- gsub(" id=\"toc-[^\"]*\"", "", parts[1])
    css <- if (length(parts) >= 3L) parts[3] else ""
    list(toc = trimws(toc), body = parts[2], css = css)
}

## ---------------------------------------------------------------------------
## Build page bodies
## ---------------------------------------------------------------------------
highlight_css <- ""
built <- list()

for (p in pages) {
    if (!is.null(p$readme)) {
        md <- character()
        for (s in p$readme) {
            body <- readme_section(s)
            if (isTRUE(p$keep_h2)) {
                md <- c(md, paste("##", s), body)
            } else {
                md <- c(md, promote(body))
            }
        }
        res <- pandoc(fix_paths(md))
    } else if (!is.null(p$file)) {
        md <- readLines(file.path(src, "pages", p$file), warn = FALSE)
        md <- gsub("@VERSION@", version, md, fixed = TRUE)
        res <- pandoc(md)
    } else if (isTRUE(p$news)) {
        md <- readLines(file.path(root, "NEWS.md"), warn = FALSE)
        fenced <- in_fence(md)
        hit <- grepl("^#{1,5} ", md) & !fenced
        md[hit] <- paste0("#", md[hit])
        res <- pandoc(md)
    } else {
        next
    }
    if (nzchar(trimws(res$css))) highlight_css <- res$css
    built[[p$id]] <- list(id = p$id, title = p$title, toc = res$toc,
                          body = res$body)
}

## ---------------------------------------------------------------------------
## Function reference from man/*.Rd
## ---------------------------------------------------------------------------
rd_files <- sort(list.files(file.path(root, "man"), pattern = "\\.Rd$",
                            full.names = TRUE))
ref <- list()
for (f in rd_files) {
    topic <- sub("\\.Rd$", "", basename(f))
    tmp <- tempfile(fileext = ".html")
    ## Rd2HTML reports help-link redirects it cannot resolve outside an
    ## installed library; they do not affect the page, so keep them quiet.
    suppressWarnings(utils::capture.output(tools::Rd2HTML(
        tools::parse_Rd(f), out = tmp, package = c("ModuleViz", version),
        no_links = TRUE), type = "message"))
    html <- paste(trimws(readLines(tmp, warn = FALSE), which = "right"),
                  collapse = "\n")
    html <- sub("(?s)^.*?<body[^>]*>", "", html, perl = TRUE)
    html <- sub("(?s)</body>.*$", "", html, perl = TRUE)
    html <- sub("(?s)<table[^>]*>.*?</table>", "", html, perl = TRUE)
    html <- gsub("</?(main|div)[^>]*>", "", html)
    title <- strip_tags(sub("(?s)^.*?<h2[^>]*>(.*?)</h2>.*$", "\\1", html,
                            perl = TRUE))
    html <- sub("(?s)<h2[^>]*>.*?</h2>", "", html, perl = TRUE)
    html <- sub("(?s)<hr ?/?>\\s*(<(div|p)[^>]*>)?\\s*\\[Package.*$", "", html,
                perl = TRUE)
    html <- sub("<hr ?/?>\\s*$", "", html, perl = TRUE)
    ## Section headings: h3 -> h2 with ids, collected for the page contents.
    heads <- regmatches(html, gregexpr("<h3[^>]*>.*?</h3>", html,
                                       perl = TRUE))[[1]]
    toc <- character()
    for (h in heads) {
        label <- strip_tags(h)
        id <- gsub("[^a-z0-9]+", "-", tolower(label))
        html <- sub(h, sprintf("<h2 id=\"%s\">%s</h2>", id, esc(label)), html,
                    fixed = TRUE)
        toc <- c(toc, sprintf("<li><a href=\"#%s\">%s</a></li>", id,
                              esc(label)))
    }
    html <- gsub("<pre><code class='language-R'>",
                 "<pre class=\"sourceCode r\"><code>", html, fixed = TRUE)
    rd_txt <- readLines(f, warn = FALSE)
    aliases <- sub("^\\\\alias\\{(.*)\\}$", "\\1",
                   grep("^\\\\alias\\{", rd_txt, value = TRUE))
    aliases <- aliases[!grepl("^print\\.", aliases)]
    ref[[topic]] <- list(
        id = paste0("reference/", topic), topic = topic, title = topic,
        heading = title, aliases = aliases,
        toc = if (length(toc)) paste0("<ul>", paste(toc, collapse = ""),
                                      "</ul>") else "",
        body = paste0("<p class=\"ref-title\">", esc(title), "</p>\n", html))
}

## Reference index: one row per help page.
rows <- vapply(ref, function(r) {
    fns <- paste(sprintf("<code>%s</code>", esc(r$aliases)), collapse = ", ")
    sprintf("<tr><td><a href=\"%s.html\">%s</a></td><td>%s</td><td>%s</td></tr>",
            r$topic, esc(r$topic), fns, esc(r$heading))
}, character(1))
built[["reference/index"]] <- list(
    id = "reference/index", title = "Function reference", toc = "",
    body = paste0(
        "<p>One page per help topic, generated from the package's ",
        "<code>man/</code> files (version ", esc(version), "). Each page ",
        "lists usage, arguments, return value, and runnable examples. The ",
        "same pages are available in R with <code>?topic</code>.</p>\n",
        "<p>For worked examples of each function, see the ",
        "<a href=\"../functions.html\">Function guide</a>.</p>\n",
        "<table><thead><tr><th>Help page</th><th>Covers</th>",
        "<th>Title</th></tr></thead><tbody>\n",
        paste(rows, collapse = "\n"), "\n</tbody></table>"))

## Final page order: reference pages follow the reference index.
order_ids <- unlist(lapply(pages, function(p) {
    if (isTRUE(p$reference_index))
        c(p$id, vapply(ref, function(r) r$id, character(1)))
    else p$id
}), use.names = FALSE)
all_pages <- c(built, stats::setNames(ref, vapply(ref, function(r) r$id,
                                                  character(1))))
all_pages <- all_pages[order_ids]

## ---------------------------------------------------------------------------
## Page shell
## ---------------------------------------------------------------------------
nav_html <- function(current, prefix) {
    in_ref <- startsWith(current, "reference/")
    items <- character()
    for (p in pages) {
        if (p$id == "index") next
        is_current <- p$id == current
        open_ref <- isTRUE(p$reference_index) && in_ref
        cls <- if (is_current) " class=\"current\"" else
            if (open_ref) " class=\"open\"" else ""
        li <- sprintf("<li%s><a href=\"%s%s.html\">%s</a>", cls, prefix, p$id,
                      esc(p$title))
        if (is_current && nzchar(all_pages[[current]]$toc))
            li <- paste0(li, all_pages[[current]]$toc)
        if (open_ref) {
            sub_items <- vapply(ref, function(r) {
                here <- r$id == current
                s <- sprintf("<li%s><a href=\"%s%s.html\">%s</a>",
                             if (here) " class=\"current\"" else "",
                             prefix, r$id, esc(r$title))
                if (here && nzchar(r$toc)) s <- paste0(s, r$toc)
                paste0(s, "</li>")
            }, character(1))
            li <- paste0(li, "<ul>", paste(sub_items, collapse = ""), "</ul>")
        }
        items <- c(items, paste0(li, "</li>"))
    }
    paste0("<ul class=\"nav\">", paste(items, collapse = "\n"), "</ul>")
}

shell <- paste(readLines(file.path(src, "page.template"), warn = FALSE),
               collapse = "\n")
fill <- function(tpl, vals) {
    for (k in names(vals))
        tpl <- gsub(paste0("{{", k, "}}"), vals[[k]], tpl, fixed = TRUE)
    tpl
}

ids <- names(all_pages)
search <- character()
for (i in seq_along(ids)) {
    p <- all_pages[[i]]
    prefix <- if (startsWith(p$id, "reference/")) "../" else ""
    is_home <- p$id == "index"
    crumbs <- if (is_home) "" else if (startsWith(p$id, "reference/") &&
                                       p$id != "reference/index")
        sprintf("<a href=\"index.html\">Function reference</a> &raquo; %s",
                esc(p$title))
    else esc(p$title)
    link <- function(j, label, cls) {
        if (j < 1 || j > length(ids)) return("")
        sprintf("<a class=\"btn %s\" href=\"%s%s.html\" title=\"%s\">%s</a>",
                cls, prefix, ids[j], esc(all_pages[[j]]$title), label)
    }
    html <- fill(shell, list(
        TITLE = if (is_home) "ModuleViz documentation" else
            paste0(esc(p$title), " &mdash; ModuleViz ", esc(version),
                   " documentation"),
        H1 = if (is_home) "" else paste0("<h1>", esc(p$title), "</h1>"),
        ROOT = prefix,
        VERSION = esc(version),
        REPO = repo,
        NAV = nav_html(p$id, prefix),
        CRUMBS = if (is_home) "" else paste0(" &raquo; ", crumbs),
        BODY = p$body,
        PREV = link(i - 1, "&larr; Previous", "prev"),
        NEXT = link(i + 1, "Next &rarr;", "next"),
        HIGHLIGHT = highlight_css,
        YEAR = format(Sys.Date(), "%Y")
    ))
    dest <- file.path(out, paste0(p$id, ".html"))
    writeLines(html, dest, useBytes = TRUE)
    search <- c(search, sprintf("{u:%s,t:%s,x:%s}",
                                js_string(paste0(p$id, ".html")),
                                js_string(p$title),
                                js_string(strip_tags(p$body))))
}
writeLines(paste0("window.MV_SEARCH=[\n", paste(search, collapse = ",\n"),
                  "\n];"),
           file.path(out, "search-index.js"), useBytes = TRUE)
invisible(file.copy(file.path(src, "site.css"), file.path(out, "site.css"),
          overwrite = TRUE))
invisible(file.copy(file.path(src, "site.js"), file.path(out, "site.js"),
          overwrite = TRUE))
## Tell GitHub Pages to serve the folder as-is (no Jekyll processing).
invisible(file.create(file.path(out, ".nojekyll")))
if (!nzchar(trimws(highlight_css))) message("note: pandoc returned no highlighting CSS")
message("built ", length(ids), " pages in ", out)
