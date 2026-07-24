#' Example longitudinal expression dataset.
#'
#' A small synthetic feature-by-sample matrix with matched sample metadata,
#' designed for examples, tests, and the package vignette. The data include 60
#' features measured across two treatment groups, four time points, and two
#' replicates per group/time combination. Features are simulated to represent
#' early-high, middle-high, late-high, and mostly flat temporal patterns.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{mat}{Numeric matrix with 60 features in rows and 16 samples in columns.}
#'   \item{meta}{Data frame with columns `SampleID`, `Condition`, `TimeHr`, and `Replicate`.}
#' }
#' @source Simulated example data generated for package tests and tutorials.
"heatmap_example"
