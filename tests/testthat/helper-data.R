example_data <- function() {
  env <- new.env()
  utils::data("heatmap_example", package = "ModuleViz", envir = env)
  env$heatmap_example
}

example_obj <- function(...) {
  ex <- example_data()
  longitudinal_cluster(ex$mat, ex$meta, sample_col = "SampleID",
                       time_col = "TimeHr", group_col = "Condition",
                       k = 4, stability = FALSE, ...)
}

## draw to a throwaway device so tests never leave Rplots.pdf behind
with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
