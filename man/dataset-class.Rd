% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\docType{class}
\name{dataset-class}
\alias{dataset-class}
\title{Training set and test set}
\description{
This class contains a dataset divided into four parts: the training set and test set, description and class labels.
}
\section{Slots}{

\describe{
\item{\code{train.x}}{the training set (description), as a \code{data.frame} or a \code{matrix}.}

\item{\code{train.y}}{the training set (target), as a \code{vector} or a \code{factor}.}

\item{\code{test.x}}{the training set (description), as a \code{data.frame} or a \code{matrix}.}

\item{\code{test.y}}{the training set (target), as a \code{vector} or a \code{factor}.}
}}

\seealso{
\code{\link{splitdata}}
}
