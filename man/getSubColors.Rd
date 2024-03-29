% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{getSubColors}
\alias{getSubColors}
\title{Run single-cell mapping}
\usage{
getSubColors(
  tab,
  class = "annotation_level_2",
  subclass = "annotation_level_4",
  external_pal = NULL,
  muted = NULL,
  pal = "Set3",
  max_pal = 12,
  random = T,
  seed = 200,
  into = "#EDEDED",
  add_n = 1
)

getSubColors(
  tab,
  class = "annotation_level_2",
  subclass = "annotation_level_4",
  external_pal = NULL,
  muted = NULL,
  pal = "Set3",
  max_pal = 12,
  random = T,
  seed = 200,
  into = "#EDEDED",
  add_n = 1
)
}
\arguments{
\item{tab}{Input data of the class data.frame}

\item{class}{Character value; The col containing the main class}

\item{subclass}{Character value; The col containing the subclass}

\item{pal}{Character value; Color pal from brewer.pal()}

\item{random}{Logical. If TRUE random mixing colors from pal}

\item{seed}{Integer value. set.seed() for constant color alignment}

\item{into}{Character value; The color for non-classified samples}

\item{add_n}{Integer value. Adopting the contrast of the subclass colors}
}
\value{
data.frame with intensity per spot
}
\description{
Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset

Get summary of image intensity per spot
}
\author{
Dieter Henrik Heiland
}
