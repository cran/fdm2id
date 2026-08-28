## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("MASS")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (hills, package = "MASS")
summary (hills)

## ----fig.height = 5-----------------------------------------------------------
plotdata (hills)

## -----------------------------------------------------------------------------
# Reproducible without a seed: leave-one-out builds n folds of one observation each, so there
# is nothing to draw. It is the one protocol of the package that needs no 'seed'.
performance (LINREG, hills [, 1], hills [, 3], protocol = "loocv", eval = c ("adjr2", "msep"))
performance (LINREG, hills [, 2], hills [, 3], protocol = "loocv", eval = c ("adjr2", "msep"))

## -----------------------------------------------------------------------------
round (cor (hills), 3)

## -----------------------------------------------------------------------------
performance (LINREG, hills [, -3], hills [, 3], protocol = "loocv", eval = c ("adjr2", "msep"))

## ----fig.height = 4.5---------------------------------------------------------
model = LINREG (hills [, -3], hills [, 3])
resplot (model)
resplot (model, index = 0)
resplot (model, index = 1)
resplot (model, index = 2)

## -----------------------------------------------------------------------------
head (sort (abs (residuals (model$model)), decreasing = TRUE), 3)
hills [c ("Knock Hill", "Bens of Jura"), ]

## -----------------------------------------------------------------------------
hills2 = cbind (hills, hills$climb^2)
colnames (hills2) = c (colnames (hills), "climb2")
performance (LINREG, hills2 [, -3], hills2 [, 3], protocol = "loocv",
             eval = c ("adjr2", "msep"))

