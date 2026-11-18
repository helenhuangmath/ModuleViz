The package includes three Bioconductor-style vignettes. Each is rendered here
with its code and figures.

| Vignette | Content |
|---|---|
| [Complete user guide](articles/ModuleViz.html) | Every function, with options and figures. |
| [Published CD8 T cell exhaustion workflow](articles/cd8-exhaustion-published-data.html) | RNA-seq modules with H3K27ac and H3K27me3 CUT&RUN aligned by gene. |
| [Longitudinal heatmap tutorial](articles/longitudinal-heatmap-tutorial.html) | Bundled and public time-course data, single-cell and multi-omics examples. |

After installation, the same vignettes open from R:

```r
vignette("ModuleViz")
vignette("cd8-exhaustion-published-data")
vignette("longitudinal-heatmap-tutorial")
```

Runnable scripts for each workflow are in the
[`examples/`](https://github.com/helenhuangmath/ModuleViz/tree/main/examples)
folder of the repository.
