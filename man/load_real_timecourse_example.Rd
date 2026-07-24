\name{load_real_timecourse_example}
\alias{load_real_timecourse_example}
\alias{available_real_timecourse_datasets}
\title{Load supported public time-course example datasets}
\description{
List and load public time-course datasets from Bioconductor
\code{timecoursedata}, returning plain matrices and metadata ready for
ModuleViz workflows.
}
\usage{
available_real_timecourse_datasets()

load_real_timecourse_example(dataset = c("sorghum_leaf", "sorghum_root",
  "influenza_mouse"), top_n_features = 1000, genotype = "BT642",
  condition = NULL, group = NULL, min_week = 3, log_transform = TRUE)
}
\arguments{
\item{dataset}{Dataset key.}
\item{top_n_features}{Number of most-variable features to keep.}
\item{genotype, condition}{Optional filters for sorghum datasets.}
\item{group}{Optional filter for influenza groups.}
\item{min_week}{Minimum sorghum week retained.}
\item{log_transform}{Use \code{log2(x + 1)} before feature selection.}
}
\value{
\code{available_real_timecourse_datasets()} returns a data table. 
\code{load_real_timecourse_example()} returns a list with \code{mat}, \code{meta},
column-name fields, source text, and a citation hint.
}
