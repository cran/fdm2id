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
data (autompg)
autompg = autompg [, -7]
summary (autompg)

## -----------------------------------------------------------------------------
apply (autompg [, -7], 2, sd)

## -----------------------------------------------------------------------------
autompg [, -7] = scale (autompg [, -7])

## ----fig.height = 6-----------------------------------------------------------
plotdata (autompg, type = "pairs", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "boxplot", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "parallel", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "histogram", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "pca", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "cda", labels = FALSE)

## ----fig.height = 5-----------------------------------------------------------
plotdata (autompg, type = "svd", labels = FALSE)

## ----fig.height = 5, eval = available && requireNamespace ("Rtsne", quietly = TRUE)----
# Variable, and visibly so: t-SNE starts from a random embedding and optimises it. Two runs
# without 'seed' give two different pictures -- same groups, different positions and shapes.
plotdata (autompg, type = "tsne", labels = FALSE, perplexity = 50, seed = 0)

