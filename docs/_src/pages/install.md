`ModuleViz` @VERSION@ is installed from GitHub. It needs R 4.3 or later.

## Install from GitHub

`BiocManager` installs the package together with its CRAN and Bioconductor
dependencies. Installing from a GitHub repository also needs the `remotes`
package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

BiocManager::install("helenhuangmath/ModuleViz")
```

Check the installation:

```r
library(ModuleViz)
packageVersion("ModuleViz")
```

## Dependencies

These are installed automatically:

| Source | Packages |
|---|---|
| Bioconductor | `ComplexHeatmap`, `SummarizedExperiment`, `S4Vectors` |
| CRAN | `circlize`, `cluster`, `data.table`, `ggplot2` (3.4.0 or later), `Matrix`, `RColorBrewer`, `withr` |

## Optional packages

Optional packages are used for single-cell input, vignettes, examples, and
package checks:

```r
BiocManager::install(c("SingleCellExperiment", "BiocStyle", "timecoursedata"))
install.packages(c("SeuratObject", "knitr", "rmarkdown", "testthat"))
```

| Package | Needed for |
|---|---|
| `SingleCellExperiment` | passing a `SingleCellExperiment` to `pseudobulk_cells()` |
| `SeuratObject` | passing a `Seurat` object to `pseudobulk_cells()` |
| `timecoursedata` | `load_real_timecourse_example()` and the public-data examples |
| `BiocStyle`, `knitr`, `rmarkdown` | building the vignettes |
| `testthat` | running the test suite |

## Development setup

To work from a checked-out copy of the repository:

```bash
git clone https://github.com/helenhuangmath/ModuleViz.git
cd ModuleViz
```

```r
devtools::load_all(".")
```

To build and check the package:

```bash
R CMD build .
R CMD check --no-manual ModuleViz_*.tar.gz
```

## Bioconductor

`ModuleViz` is not on Bioconductor yet. Until it is,
`BiocManager::install("ModuleViz")` will not find it; use the GitHub form
above.
