% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{intern.intraclass}
\alias{intern.intraclass}
\title{Clustering evaluation through intraclass inertia}
\usage{
intern.intraclass(clus, d, type = c("global", "cluster"))
}
\arguments{
\item{clus}{The extracted clusters.}

\item{d}{The dataset.}

\item{type}{Indicates whether a "global" or a "cluster"-wise evaluation should be used.}
}
\value{
The evaluation of the clustering.
}
\description{
Evaluation a clustering algorithm according to intraclass inertia.
}
\examples{
require (datasets)
data (iris)
km = KMEANS (iris [, -5], k = 3)
intern.intraclass (km$clus, iris [, -5])
}
\seealso{
\code{\link{intern}}, \code{\link{intern.dunn}}, \code{\link{intern.interclass}}
}
