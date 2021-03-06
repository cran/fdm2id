% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\name{data.xor}
\alias{data.xor}
\title{XOR dataset}
\usage{
data.xor(
  n = 100,
  ndim = 2,
  sigma = 0.25,
  levels = NULL,
  graph = TRUE,
  seed = NULL
)
}
\arguments{
\item{n}{Number of observations in each cluster.}

\item{ndim}{The number of dimensions (2^ndim clusters are formed, grouped into two classes).}

\item{sigma}{The variance.}

\item{levels}{Name of each class.}

\item{graph}{A logical indicating whether or not a graphic should be plotted.}

\item{seed}{A specified seed for random number generation.}
}
\value{
A randomly generated dataset.
}
\description{
Generate "XOR" dataset.
}
\examples{
data.xor ()
}
\seealso{
\code{\link{data.diag}}, \code{\link{data.gauss}}, \code{\link{data.parabol}}, \code{\link{data.target2}}, \code{\link{data.twomoons}}
}
