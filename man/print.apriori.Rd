% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/association.R
\name{print.apriori}
\alias{print.apriori}
\title{Print a classification model obtained by APRIORI}
\usage{
\method{print}{apriori}(x, ...)
}
\arguments{
\item{x}{The model to be printed.}

\item{...}{Other parameters.}
}
\description{
Print the set of rules in the classification model.
}
\examples{
require ("datasets")
data (iris)
d = discretizeDF (iris,
    default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
print (model)
}
\seealso{
\code{\link{APRIORI}}, \code{\link{predict.apriori}}, \code{\link{summary.apriori}},
\code{\link{apriori-class}}, \code{\link[arules]{apriori}}
}
