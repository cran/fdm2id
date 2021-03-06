% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{MLP}
\alias{MLP}
\title{Classification using Multilayer Perceptron}
\usage{
MLP(
  train,
  labels,
  hidden = ifelse(is.vector(train), 2:(1 + nlevels(labels)), 2:(ncol(train) +
    nlevels(labels))),
  decay = 10^(-3:-1),
  methodparameters = NULL,
  tune = FALSE,
  ...
)
}
\arguments{
\item{train}{The training set (description), as a \code{data.frame}.}

\item{labels}{Class labels of the training set (\code{vector} or \code{factor}).}

\item{hidden}{The size of the hidden layer (if a vector, cross-over validation is used to chose the best size).}

\item{decay}{The decay (between 0 and 1) of the backpropagation algorithm (if a vector, cross-over validation is used to chose the best size).}

\item{methodparameters}{Object containing the parameters. If given, it replaces \code{size} and \code{decay}.}

\item{tune}{If true, the function returns paramters instead of a classification model.}

\item{...}{Other parameters.}
}
\value{
The classification model.
}
\description{
This function builds a classification model using Multilayer Perceptron.
}
\examples{
\dontrun{
require (datasets)
data (iris)
MLP (iris [, -5], iris [, 5], hidden = 4, decay = .1)
}
}
\seealso{
\code{\link[nnet]{nnet}}
}
