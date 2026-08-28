## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("mlbench")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
library (mlbench)
data (BreastCancer)
BreastCancer = BreastCancer [, -1]
BreastCancer [, -10] = lapply (BreastCancer [, -10], function (x) as.numeric (as.character (x)))
names (which (colSums (is.na (BreastCancer)) > 0))
train = BreastCancer [!is.na (BreastCancer$Bare.nuclei), ]
BreastCancer [is.na (BreastCancer$Bare.nuclei), 6] =
  round (predict (LINREG (train [, -c (6, 10)], train [, 6]),
                  BreastCancer [is.na (BreastCancer$Bare.nuclei), -6]))
summary (BreastCancer)

## ----fig.height = 6-----------------------------------------------------------
plotdata (BreastCancer)

## -----------------------------------------------------------------------------
pca = PCA (BreastCancer, quali.sup = 10, scale.unit = TRUE)
kaiser (pca)

## -----------------------------------------------------------------------------
plot (pca, type = "eig")

## ----fig.height = 5.5---------------------------------------------------------
plot (pca, type = "cor")

## -----------------------------------------------------------------------------
round (pca$var$coord [, 1:2], 2)

## ----fig.height = 5-----------------------------------------------------------
plotdata (pca$ind$coord [, 1:2], BreastCancer [, 10], type = "scatter")

## -----------------------------------------------------------------------------
round (tapply (pca$ind$coord [, 1], BreastCancer [, 10], mean), 2)
round (tapply (pca$ind$coord [, 1], BreastCancer [, 10], sd), 2)

