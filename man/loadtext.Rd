% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text.R
\name{loadtext}
\alias{loadtext}
\title{load a text file}
\usage{
loadtext(
  file = file.choose(),
  dir = "~/",
  collapse = TRUE,
  sep = NULL,
  categories = NULL
)
}
\arguments{
\item{file}{The path or URL of the text file.}

\item{dir}{The (temporary) directory, where the file is downloaded. The file is deleted at the end of this function.}

\item{collapse}{Indicates whether or not lines of each documents should collapse together or not.}

\item{sep}{Separator between text fields.}

\item{categories}{Columns that should be considered as categorial data.}
}
\value{
The text contained in the dowloaded file.
}
\description{
(Down)Load a text file (and extract it if it is in a zip file).
}
\examples{
\dontrun{
text = loadtext ("http://mattmahoney.net/dc/text8.zip")
}
}
\seealso{
\code{\link[utils]{download.file}}, \code{\link[utils]{unzip}}
}
