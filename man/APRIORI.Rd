% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/association.R
\name{APRIORI}
\alias{APRIORI}
\title{Classification using APRIORI}
\usage{
APRIORI(
  train,
  labels,
  supp = 0.05,
  conf = 0.8,
  prune = FALSE,
  tune = FALSE,
  ...
)
}
\arguments{
\item{train}{The training set (description), as a \code{data.frame}.}

\item{labels}{Class labels of the training set (\code{vector} or \code{factor}).}

\item{supp}{The minimal support of an item set (numeric value).}

\item{conf}{The minimal confidence of an item set (numeric value).}

\item{prune}{A logical indicating whether to prune redundant rules or not (default: \code{FALSE}).}

\item{tune}{If true, the function returns paramters instead of a classification model.}

\item{...}{Other parameters.}
}
\value{
The classification model, as an object of class \code{apriori}.
}
\description{
This function builds a classification model using the association rules method APRIORI.
}
\examples{
require ("datasets")
data (iris)
d = discretizeDF (iris,
    default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
}
\seealso{
\code{\link{predict.apriori}}, \code{\link{apriori-class}}, \code{\link[arules]{apriori}}
}
