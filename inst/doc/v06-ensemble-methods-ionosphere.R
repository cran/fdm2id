## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("glmnet", "rpart", "rpart.plot", "randomForest", "xgboost", "ROCR")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (ionosphere)
dim (ionosphere)
table (ionosphere [, 34])

## ----fig.height = 5-----------------------------------------------------------
plotdata (ionosphere [, -34], ionosphere [, 34], type = "pca")

## -----------------------------------------------------------------------------
# Variable on both counts: the bootstrap draws the ten resamples, and four of the six methods
# are randomised in themselves -- bagging and the forest draw their samples and their
# variables, the two boosting methods their subsamples. Without 'seed' the table moves by a
# few thousandths at every run, which is more than the gap between the two leaders.
performance (c (LR, CART, BAGGING, RANDOMFOREST, ADABOOST, GRADIENTBOOSTING),
             ionosphere [, -34], ionosphere [, 34], type = "evaluation",
             protocol = "bootstrap", eval = "accuracy", nruns = 10, seed = 0,
             learningmethod = CART)

## ----fig.height = 5-----------------------------------------------------------
performance (c (RANDOMFOREST, ADABOOST), ionosphere [, -34], ionosphere [, 34],
             type = "roc", protocol = "bootstrap", nruns = 10, fuzzy = TRUE, seed = 0,
             learningmethod = CART)

## -----------------------------------------------------------------------------
performance (c (RANDOMFOREST, ADABOOST), ionosphere [, -34], ionosphere [, 34],
             type = "evaluation", protocol = "bootstrap",
             eval = c ("precision", "recall"), nruns = 10, seed = 0, learningmethod = CART)

## ----fig.height = 4.5---------------------------------------------------------
performance (RANDOMFOREST, ionosphere [, -34], ionosphere [, 34], type = "confusion",
             protocol = "bootstrap", nruns = 10, seed = 0)
performance (ADABOOST, ionosphere [, -34], ionosphere [, 34], type = "confusion",
             protocol = "bootstrap", nruns = 10, seed = 0, learningmethod = CART)

