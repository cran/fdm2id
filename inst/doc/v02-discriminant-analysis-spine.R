## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("MASS", "e1071", "glmnet", "rpart", "rpart.plot")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (spine)
summary (spine)

## ----fig.height = 6-----------------------------------------------------------
plotdata (spine, k = spine [, 7])

## ----fig.height = 6-----------------------------------------------------------
plotdata (spine, k = spine [, 8])

## -----------------------------------------------------------------------------
# Variable: the bootstrap draws its 100 resamples at random, so without 'seed' this table
# changes at every run -- by a few thousandths here, enough to swap two close methods.
performance (c (NB, LDA, CDA, LR), spine [, 1:6], spine [, 7], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
performance (c (NB, LDA, CDA, LR), spine [, 1:6], spine [, 8], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
performance (LR, spine [, 1:6], spine [, 7], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)
performance (LR, spine [, 1:6], spine [, 8], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)

## ----fig.height = 4.5---------------------------------------------------------
performance (LR, spine [, 1:6], spine [, 8], type = "confusion",
             protocol = "bootstrap", nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
performance (c (NB, LDA, CDA, LR, SVMl), spine [, 1:6], spine [, 8], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)

## -----------------------------------------------------------------------------
# Variable, twice over: on top of the bootstrap, KNN, MLP and the SVMs search their
# hyperparameter grid by cross-validation, so the model itself is drawn at random too.
performance (c (KNN, CART, MLP, SVMr), spine [, 1:6], spine [, 8], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 100, seed = 0)

## ----fig.height = 5-----------------------------------------------------------
cartplot (CART (spine [, 1:6], spine [, 8]))

