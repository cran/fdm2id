## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set (collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.2,
                       fig.align = "center")
optional = c ("cluster")
available = all (sapply (optional, requireNamespace, quietly = TRUE))
knitr::opts_chunk$set (eval = available)

## ----echo = FALSE, eval = !available, results = "asis"------------------------
# cat ("**Note.** This vignette needs the following packages, some of which are missing:",
#      paste (optional, collapse = ", "), "-- the code is shown but not run.\n")

## ----message = FALSE, warning = FALSE-----------------------------------------
library (fdm2id)

## -----------------------------------------------------------------------------
data (wheat)
summary (wheat [, -8])

## ----fig.height = 6-----------------------------------------------------------
plotdata (wheat [, -8])

## -----------------------------------------------------------------------------
apply (wheat [, -8], 2, sd)

## -----------------------------------------------------------------------------
wheat [, -8] = scale (wheat [, -8])

## -----------------------------------------------------------------------------
# Variable: K-means starts from centres drawn at random. 'nstart = 100' keeps the best of a
# hundred starts, which makes the answer stable in practice, but only 'seed' makes it exact.
kmeans.getk (wheat [, -8], nstart = 100, graph = TRUE, seed = 0)

## ----fig.height = 5-----------------------------------------------------------
single = HCA (wheat [, -8], method = "single")
plotclus (single, wheat, "tree")

## ----fig.height = 5-----------------------------------------------------------
ward = HCA (wheat [, -8], method = "ward")
plotclus (ward, wheat, "tree")

## -----------------------------------------------------------------------------
sapply (c ("pseudo-F", "silhouette", "elbow"),
        function (criterion) kmeans.getk (wheat [, -8], criterion = criterion,
                                          nstart = 100, seed = 0))

## -----------------------------------------------------------------------------
km = KMEANS (wheat [, -8], k = 3, nstart = 100, seed = 0)
ward = HCA (wheat [, -8], k = 3, method = "ward")
intern (km, wheat [, -8], eval = c ("intraclass", "interclass"))
intern (ward, wheat [, -8], eval = c ("intraclass", "interclass"))

## -----------------------------------------------------------------------------
# Variable: stability resamples the dataset. Note that HCA itself is deterministic -- here the
# randomness is entirely in the resampling, not in the method being judged.
stability (KMEANS, wheat [, -8], type = "global", k = 3, nstart = 100, seed = 0)
stability (HCA, wheat [, -8], type = "global", method = "ward", k = 3, seed = 0)

## -----------------------------------------------------------------------------
compare (km, wheat [, 8], comp = "pairwise")
compare (ward, wheat [, 8], comp = "pairwise")

## -----------------------------------------------------------------------------
table (km$cluster, wheat [, 8])
table (ward$cluster, wheat [, 8])

