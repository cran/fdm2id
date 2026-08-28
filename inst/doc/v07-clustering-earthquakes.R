## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("flexclust", "RSpectra")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (quakes, package = "datasets")
summary (quakes [, -5])

## ----fig.height = 6-----------------------------------------------------------
plotdata (quakes [, -5])

## ----fig.height = 5-----------------------------------------------------------
plotdata (quakes [, 1:2], type = "scatter")

## -----------------------------------------------------------------------------
# Variable, all three: K-means starts from random centres, EM from a random initialisation,
# and spectral clustering ends on a K-means in the eigenvector space. Without 'seed', the
# cluster numbering alone changes from run to run -- so every "cluster 3" below would have to
# be re-read against the picture.
km = KMEANS (quakes [, 1:2], k = 3, seed = 0)
em = EM (quakes [, 1:2], k = 3, seed = 0)
sc = SPECTRAL (quakes [, 1:2], k = 3, sigma = .5, seed = 0)

## -----------------------------------------------------------------------------
# Variable: the resampling, plus the method re-run on each resample. The same 'seed' is given
# to the three calls so that they compare the methods on the same resamples.
stability (KMEANS, quakes [, 1:2], km, k = 3, seed = 0)
stability (EM, quakes [, 1:2], em, k = 3, seed = 0)
stability (SPECTRAL, quakes [, 1:2], sc, k = 3, sigma = .5, seed = 0)

## ----fig.height = 4-----------------------------------------------------------
plotdata (quakes [, "depth"], km$cluster, type = "boxplot")
plotdata (quakes [, "depth"], em$cluster, type = "boxplot")
plotdata (quakes [, "depth"], sc$cluster, type = "boxplot")

## ----fig.height = 4-----------------------------------------------------------
plotdata (quakes [, "mag"], km$cluster, type = "boxplot")
plotdata (quakes [, "mag"], em$cluster, type = "boxplot")
plotdata (quakes [, "mag"], sc$cluster, type = "boxplot")

## ----fig.height = 5-----------------------------------------------------------
plotclus (km, quakes [, 1:2])
plotclus (em, quakes [, 1:2])
plotclus (sc, quakes [, 1:2])

## ----fig.height = 5-----------------------------------------------------------
plotdata (quakes [, 2:1], cut (quakes$depth, c (0, 150, 350, 700)), type = "scatter")

## -----------------------------------------------------------------------------
west = quakes [quakes$long < 176, ]
table (deep = west$depth > 350)
round (apply (west [west$depth > 350, 1:2], 2, range), 1)

