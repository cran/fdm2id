## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("MASS", "e1071")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (crabs, package = "MASS")
summary (crabs)

## ----fig.height = 6-----------------------------------------------------------
plotdata (crabs [, 4:8], crabs [, 1], type = "pairs")

## ----fig.height = 6-----------------------------------------------------------
plotdata (crabs [, 4:8], crabs [, 2], type = "pairs")

## -----------------------------------------------------------------------------
round (cor (crabs [, 4:8]), 3)

## -----------------------------------------------------------------------------
# The two filters are deterministic; the wrapper is not -- it judges each subset by fitting a
# naive Bayes classifier under a bootstrap, so without 'seed' it can stop at a different
# subset from one run to the next. That is the criterion's own variance, not the data's.
s.fstat1 = selectfeatures (crabs [, 4:8], crabs [, 1], algorithm = "ranking",
                           unieval = "fisher", multieval = "fstat", seed = 0)
s.mrmr1 = selectfeatures (crabs [, 4:8], crabs [, 1], algorithm = "ranking",
                          unieval = "fisher", multieval = "mrmr", seed = 0)
s.wrap1 = selectfeatures (crabs [, 4:8], crabs [, 1], algorithm = "ranking",
                          unieval = "fisher", multieval = "wrapper", wrapmethod = NB,
                          seed = 0)
s.fstat1
s.mrmr1
s.wrap1

## -----------------------------------------------------------------------------
performance (NB, crabs [, 4:8], crabs [, 1], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.fstat1$selection], crabs [, 1], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.mrmr1$selection], crabs [, 1], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.wrap1$selection], crabs [, 1], nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
s.fstat2 = selectfeatures (crabs [, 4:8], crabs [, 2], algorithm = "ranking",
                           unieval = "fisher", multieval = "fstat", seed = 0)
s.mrmr2 = selectfeatures (crabs [, 4:8], crabs [, 2], algorithm = "ranking",
                          unieval = "fisher", multieval = "mrmr", seed = 0)
s.wrap2 = selectfeatures (crabs [, 4:8], crabs [, 2], algorithm = "ranking",
                          unieval = "fisher", multieval = "wrapper", wrapmethod = NB,
                          seed = 0)
s.fstat2
s.mrmr2
s.wrap2

## -----------------------------------------------------------------------------
performance (NB, crabs [, 4:8], crabs [, 2], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.fstat2$selection], crabs [, 2], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.mrmr2$selection], crabs [, 2], nruns = 100, seed = 0)
performance (NB, crabs [, 4:8] [, s.wrap2$selection], crabs [, 2], nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
colnames (crabs) [4:8] [s.fstat1$selection]
colnames (crabs) [4:8] [s.wrap2$selection]

