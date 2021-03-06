% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_enveloppe.r
\name{return.mean.enveloppe}
\alias{return.mean.enveloppe}
\title{return.mean.enveloppe}
\usage{
return.mean.enveloppe(extrema, data, zcol = "z", method = "splines",
  n.pts.spline = 3, verbose = TRUE)
}
\arguments{
\item{extrema}{.}

\item{data}{.}

\item{zcol}{Name of the column containing the data.}

\item{method}{Interpolation method. Currently only `splines` is supported.}

\item{n.pts.spline}{.}

\item{verbose}{Prints progress information messages. Defaults to TRUE.}
}
\description{
Internal function, returns the mean enveloppe of a spatial dataset.
}
\author{
Pierre Roudier
}
