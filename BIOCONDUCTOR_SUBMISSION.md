# ModuleViz Bioconductor Submission Notes

This document prepares `ModuleViz` for submission to Bioconductor as a software
package. It is written for the maintainer and can be used as the basis for the
GitHub issue in `Bioconductor/Contributions`.

## Package Summary

`ModuleViz` is an R package for longitudinal and time-course omics
visualization. It clusters rows of an expression, accessibility, or other omics
matrix into temporal modules, orders modules by trajectory, and creates
publication-ready heatmaps and module pattern plots. It also supports paired
omics visualization, such as RNA and ATAC/CUT&Run data aligned by gene.

The package is intended for high-throughput genomic data analysis and
visualization, including:

- RNA-seq time courses
- ATAC-seq or CUT&Run time courses
- paired RNA/chromatin assays
- pseudobulk single-cell time-course summaries
- treatment, differentiation, or exhaustion trajectories

Core user-facing functions include:

- `longitudinal_cluster()`: cluster features into ordered temporal modules
- `longitudinal_heatmap()`: draw ordered module heatmaps using
  `ComplexHeatmap`
- `pattern_lineplot()`: summarize module trajectories as line plots
- `dual_omics_heatmap()`: align and compare two omics layers by gene
- `module_correlation_heatmap()`: inspect within-module correlation structure
- `module_correlation_summary()`: summarize within-module correlations
- `choose_k()`: help choose the number of clusters
- `write_memberships()`: export module membership, centroid, and stability
  tables
- `run_longitudinal()`: run clustering, heatmap plotting, line plotting, and
  membership export in one call

## Fit For Bioconductor

Bioconductor packages should support high-throughput genomic analysis,
interoperate with existing Bioconductor infrastructure, provide reproducible
documentation, and include evaluated examples and vignettes. `ModuleViz` fits
this scope because it provides visualization and clustering workflows for
genomic time-course assays and builds directly on Bioconductor infrastructure,
especially `ComplexHeatmap`.

The current `DESCRIPTION` already includes a `biocViews` field:

```text
Software, GeneExpression, Transcriptomics, RNASeq, TimeCourse, Clustering,
Visualization
```

Before submission, confirm these terms are valid leaf terms for the current
Bioconductor devel branch. Bioconductor requires `biocViews` to be
case-sensitive, start with lower-case `biocViews`, and include at least two
leaf-node terms from the same package type.

## Repository Requirements

Submit the default branch of the GitHub repository. Bioconductor uses the
default branch for review, so the default branch should contain only package
code and package-related files.

Recommended repository state before submission:

- Default branch contains the R package root directly.
- No generated check directories are committed, such as `ModuleViz.Rcheck` or
  `..Rcheck`.
- No local IDE folders are committed, such as `.Rproj.user`.
- No local assistant/tooling folders are committed.
- No example output directories are committed unless they are intentional
  package assets.
- Large files are avoided unless they are essential and documented.
- Vignettes are in `vignettes/` and can be evaluated.
- Tests are in `tests/testthat/`.
- The maintainer in `DESCRIPTION` is the same person submitting the package.

Current cleanup status:

- Local assistant/tooling artifacts have been removed from this workspace.
- Searches for removed-tooling names outside `.git` returned no matches.
- `.Rbuildignore` excludes common local/build artifacts:
  `.git`, `.gitignore`, `.Rproj.user`, `*.Rproj`, `results`,
  `examples/outputs`, package tarballs, `ModuleViz.Rcheck`, and `..Rcheck`.

## DESCRIPTION Checklist

Before submission, review `DESCRIPTION` for the following:

- `Package`, `Title`, `Version`, `Authors@R`, `Description`, `License`,
  `Encoding`, `Depends`, `Imports`, `Suggests`, `biocViews`, and
  `VignetteBuilder` are present and accurate.
- The package version follows Bioconductor development conventions. For a new
  package under review, use a development version such as `0.99.0`, then bump
  the version when responding to build or review changes.
- `biocViews` contains valid current Bioconductor terms.
- Runtime dependencies are in `Imports`, not `Suggests`.
- Vignette-only and example-only dependencies are in `Suggests`.
- `SystemRequirements` is present only if non-R system software is required.
- The package does not attempt to install dependencies from R code, examples,
  vignettes, or `.onLoad()`.

Current package-specific notes:

- `ComplexHeatmap` and `circlize` are core dependencies and belong in
  `Imports`.
- `BiocStyle`, `knitr`, and `rmarkdown` are appropriate vignette dependencies
  in `Suggests`.
- `timecoursedata` appears to be used for public-data examples and should remain
  in `Suggests` unless required at runtime.
- Consider changing `Version: 0.1.0` to `Version: 0.99.0` before formal
  Bioconductor submission.

## Documentation Checklist

Bioconductor review expects complete documentation and reproducible examples.

Confirm the following before submission:

- Every exported function has a matching `.Rd` help page.
- All parameters, return values, and side effects are documented.
- Examples are runnable, fast, and do not write into the user's working
  directory unless wrapped in `tempdir()` or clearly guarded.
- Vignettes build from source and use evaluated code where practical.
- Vignettes use `BiocStyle::html_document` or another Bioconductor-appropriate
  format.
- The README explains installation, core workflow, and package scope.
- The package has `NEWS.md` entries for user-facing changes.

Current documentation assets:

- `README.md`
- `NEWS.md`
- `vignettes/ModuleViz.Rmd`
- `vignettes/cd8-exhaustion-published-data.Rmd`
- `vignettes/longitudinal-heatmap-tutorial.Rmd`
- `man/*.Rd`
- runnable scripts in `examples/`

## Testing And Checks

Run these checks from a clean environment before submission:

```bash
R CMD build .
R CMD check --no-manual ModuleViz_*.tar.gz
```

If suggested packages are intentionally unavailable during a quick local check:

```bash
env _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual ModuleViz_*.tar.gz
```

Run Bioconductor-specific checks:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("BiocCheck")
BiocCheck::BiocCheck(".")
```

Also run the test suite:

```r
testthat::test_dir("tests/testthat")
```

Current local validation note:

`R CMD check --no-manual .` was attempted in this workspace, but the check could
not proceed because this local R environment does not have required packages
installed:

- `circlize`
- `ComplexHeatmap`
- `ggplot2`

It also lacks suggested packages used for a complete check:

- `BiocStyle`
- `timecoursedata`

Run the checks again after installing required and suggested dependencies in a
fresh R/Bioconductor environment.

## Bioconductor Submission Process

Submit by opening a new issue in the official
`Bioconductor/Contributions` GitHub repository. The issue should link to the
default branch of the package repository. Bioconductor review expects the
submitter to be the package maintainer listed in `DESCRIPTION`.

After initial review, Bioconductor may add the package to
`git.bioconductor.org` during the review process. Further review builds may be
triggered by version bumps and pushes to the Bioconductor git repository. Follow
reviewer instructions exactly and respond with concise summaries of changes.

Expected review cycle:

- initial submission issue is opened
- automated checks run
- a reviewer comments on technical, documentation, style, or scientific issues
- maintainer updates the package and bumps the version
- additional builds run
- accepted packages enter Bioconductor devel first
- the package becomes available to release users after the next Bioconductor
  release cycle

## Draft Submission Issue

Use the following text as a starting point for the
`Bioconductor/Contributions` issue.

```markdown
Package: ModuleViz

Repository: https://github.com/helenhuangmath/ModuleViz

Type: Software

biocViews: Software, GeneExpression, Transcriptomics, RNASeq, TimeCourse,
Clustering, Visualization

Maintainer: Hua Huang <helenhuang.math@gmail.com>

Description:
ModuleViz provides visualization and clustering tools for longitudinal
high-throughput omics data. It clusters features into ordered temporal modules,
draws ComplexHeatmap-based module heatmaps, summarizes module trajectories with
line plots, and supports paired RNA/chromatin visualizations aligned by gene.
The package is intended for RNA-seq, ATAC-seq, CUT&Run, pseudobulk single-cell,
and related genomic time-course workflows.

Why Bioconductor:
The package targets high-throughput genomic time-course analysis and depends on
Bioconductor visualization infrastructure through ComplexHeatmap. It provides
reproducible examples and vignettes for genomic use cases and is designed to
fit into existing R/Bioconductor workflows using matrix-like omics data and
sample metadata.

Major functionality:
- longitudinal_cluster(): cluster features into temporal modules
- longitudinal_heatmap(): draw ordered module heatmaps
- pattern_lineplot(): plot mean module trajectories
- dual_omics_heatmap(): compare paired omics layers by gene
- module_correlation_heatmap(): inspect within-module correlation structure
- choose_k(): generate clustering diagnostics
- write_memberships(): export module assignments and summaries

Package status:
- R package structure is complete.
- Exported functions are documented.
- Vignettes and examples are included.
- testthat tests are included.
- Local/tooling artifacts have been removed.

Additional notes:
This is a new Bioconductor software package submission. I am the package
maintainer listed in DESCRIPTION and will support the package through the
Bioconductor review process and support channels.
```

## Pre-Submission Action Items

Complete these before opening the submission issue:

- Install all required and suggested dependencies in a current Bioconductor
  devel environment.
- Run `R CMD build`.
- Run `R CMD check --no-manual` on the built tarball.
- Run `BiocCheck::BiocCheck(".")`.
- Confirm that all vignettes build successfully.
- Confirm that `examples/validate_real_data.R` is either intended for package
  distribution or move it to a non-default branch / ignore it if it is only a
  local validation script.
- Confirm that committed PNG figures are necessary package assets and not
  oversized generated output.
- Confirm `biocViews` terms against the current Bioconductor vocabulary.
- Consider changing version to `0.99.0` for new package submission.
- Merge the current PR into `main` only after review of the package changes.

## Official References

- Bioconductor package submissions:
  https://contributions.bioconductor.org/bioconductor-package-submissions.html
- Bioconductor package guidelines and review documentation:
  https://contributions.bioconductor.org/
- DESCRIPTION file guidance:
  https://contributions.bioconductor.org/description.html
- Important Bioconductor package development features:
  https://contributions.bioconductor.org/important-bioconductor-package-development-features.html
- Bioconductor Contributions issue tracker:
  https://github.com/Bioconductor/Contributions
