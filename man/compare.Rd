% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{compare}
\alias{compare}
\title{Comparison of two sets of clusters}
\usage{
compare(clus, gt, eval = "accuracy", comp = c("max", "pairwise", "cluster"))
}
\arguments{
\item{clus}{The extracted clusters.}

\item{gt}{The real clusters.}

\item{eval}{The evluation criterion.}

\item{comp}{Indicates whether a "max" or a "pairwise" evaluation should be used, or the evaluation for each individual "cluster".}
}
\value{
A numeric value indicating how much the two sets of clusters are similar.
}
\description{
Comparison of two sets of clusters
}
\examples{
require (datasets)
data (iris)
km = KMEANS (iris [, -5], k = 3)
compare (km$cluster, iris [, 5])
\dontrun{
compare (km$cluster, iris [, 5], eval = c ("accuracy", "kappa"), comp = "pairwise")
}
}
\seealso{
\code{\link{compare.accuracy}}, \code{\link{compare.jaccard}}, \code{\link{compare.kappa}}, \code{\link{intern}}, \code{\link{stability}}
}
