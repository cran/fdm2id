# Tests for clustering.R

test_that ("KMEANS / HCA smoke tests", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  expect_equal (length (km$cluster), nrow (iris))
  hca = HCA (iris [, -5], k = 3, method = "ward")
  expect_equal (length (unique (hca$cluster)), 3)
})

# --- Regression: DBSCAN's epsilonDist parameter renamed to eps (breaking change) ---------

test_that ("DBSCAN's eps parameter (renamed from epsilonDist) works", {
  skip_if_not_installed ("fpc")
  data (iris)
  model = DBSCAN (iris [, -5], minpts = 5, eps = 1)
  expect_equal (length (model$cluster), nrow (iris))
})

# --- Regression: compare()/intern() dispatch tables (get(paste(...)) -> named list) ------

test_that ("compare() is exported and its dispatch table resolves compare.* functions", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  res = compare (km$cluster, iris [, 5])
  expect_true (is.numeric (res))
  res2 = compare (km$cluster, iris [, 5], eval = c ("accuracy", "jaccard", "kappa"))
  expect_named (res2, c ("accuracy", "jaccard", "kappa"))
})

test_that ("compare() gives a clear error for an unknown criterion", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  expect_error (compare (km$cluster, iris [, 5], eval = "not_a_real_criterion"),
               "unknown evaluation criterion")
})

test_that ("intern() is exported and its dispatch table resolves intern.* functions", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  res = intern (km$cluster, iris [, -5])
  expect_true (is.numeric (res))
  res2 = intern (km$cluster, iris [, -5], eval = c ("intraclass", "interclass"))
  expect_named (res2, c ("intraclass", "interclass"))
})

test_that ("intern() gives a clear error for an unknown criterion", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  expect_error (intern (km$cluster, iris [, -5], eval = "not_a_real_criterion"),
               "unknown evaluation criterion")
})

# --- Regression: plotclus() extractor dispatch table (repeated "x" %in% method chain) ----

test_that ("plotclus() is exported and its extractor dispatch table works for kmeans and hca", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, nstart = 5)
  ward = HCA (iris [, -5], k = 3, method = "ward")
  grDevices::pdf (NULL) # discard graphics output, avoid writing Rplots.pdf during tests
  on.exit (grDevices::dev.off ())
  expect_error (plotclus (km, iris [, -5], type = "scatter"), NA)
  expect_error (plotclus (km, iris [, -5], type = "boxplot"), NA)
  expect_error (plotclus (ward, iris [, -5], type = "tree"), NA)
})

# --- Regression: accuracy0() on a constant cluster indicator -----------------------------
# accuracy0() used sum (diag (table (clus, gt))): when one of the two indicators is constant
# the table is 1 x 2 and the diagonal reads the wrong cells. Same root cause as
# eval.accuracy() in classification.R; both now go through align.labels().

test_that ("compare.accuracy() copes with a cluster that covers every observation", {
  data (iris)
  # A single cluster: every "clus == i" indicator is constant.
  expect_error (compare.accuracy (rep (1, nrow (iris)), iris [, 5]), NA)
  expect_true (is.finite (compare.accuracy (rep (1, nrow (iris)), iris [, 5])))
})

test_that ("compare() on a perfect clustering returns an accuracy of 1", {
  data (iris)
  clus = as.numeric (iris [, 5]) # the clustering *is* the ground truth
  expect_equal (unname (compare (clus, iris [, 5], eval = "accuracy")), 1)
})

# --- Regression: plot.som() on a two-variable dataset ------------------------------------
# centers.coord was initialised to x$som$codes, which is a *list*; only the "more than two
# variables" branch replaced it with codes [[1]]. And the mapping branch passed
# classif = x$unit.classif, a field that lives on the wrapped kohonen map, not on the
# fdm2id object, so it was always NULL.

test_that ("plot.som() works whatever the number of variables", {
  skip_if_not_installed ("kohonen")
  data (iris)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  som2 = SOM (iris [, 3:4], xdim = 4, ydim = 4, rlen = 100, post = "ward", k = 3, seed = 1)
  som4 = SOM (iris [, -5],  xdim = 4, ydim = 4, rlen = 100, post = "ward", k = 3, seed = 1)
  expect_error (plot (som2), NA)                    # used to fail: "incorrect number of dimensions"
  expect_error (plot (som4), NA)
  expect_error (plot (som2, type = "mapping"), NA)
  expect_error (plot (som4, type = "mapping"), NA)
})

test_that ("SOM objects carry the unit assignments plot.som() needs", {
  skip_if_not_installed ("kohonen")
  data (iris)
  som = SOM (iris [, -5], xdim = 4, ydim = 4, rlen = 100, seed = 1)
  expect_false (is.null (som$som$unit.classif))
  expect_equal (length (som$som$unit.classif), nrow (iris))
})

# --- Regression: compare() with a cluster id that is never used --------------------------
# clusters <- min (kk1):max (kk1) enumerated *possible* ids while table (kk1) counts only the
# existing ones, so weighted.mean() failed with "'x' and 'w' must have the same length".

test_that ("compare() copes with a gap in the cluster ids", {
  data (iris)
  clus = c (rep (1, 50), rep (3, 100)) # no cluster 2
  for (e in c ("accuracy", "jaccard", "kappa"))
    expect_error (compare (clus, iris [, 5], eval = e), NA, info = e)
  expect_error (compare (clus, iris [, 5], comp = "cluster"), NA)
})

test_that ("compare() is unchanged on contiguous cluster ids", {
  data (iris)
  clus = as.numeric (iris [, 5])
  expect_equal (unname (compare (clus, iris [, 5], eval = "accuracy")), 1)
})

# --- Regression: cluster centers were projected with the wrong transform -----------------
# scatterplot() centred them on the mean of the *centers* instead of the mean of the data,
# and plot.som() divided the codebook vectors by the standard deviations although the PCA
# behind the plot is run with scale.unit = FALSE. pca.project2d() now exposes the projection
# itself, so both callers use the same one.

test_that ("pca.project2d()'s projection reproduces the plotted coordinates", {
  data (iris)
  proj = fdm2id:::pca.project2d (iris [, -5])
  expect_equal (unname (proj$project (iris [, -5]) [, 1:2]), unname (proj$coord [, 1:2]))
})

test_that ("a projected cluster center is the center of its projected points", {
  data (iris)
  set.seed (1)
  km = stats::kmeans (iris [, -5], 3)
  proj = fdm2id:::pca.project2d (iris [, -5])
  centers = proj$project (km$centers)
  for (i in 1:3)
    expect_equal (unname (centers [i, 1:2]),
                  unname (colMeans (proj$coord [km$cluster == i, 1:2])))
})

test_that ("scatterplot() draws ellipses instead of failing", {
  skip_if_not_installed ("car")
  skip_if_not_installed ("mclust")
  data (iris)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  set.seed (1)
  km = stats::kmeans (iris [, -5], 3)
  # car::ellipse() used to be called with the covariance matrix as its first positional
  # argument -- i.e. as 'center' -- so this branch had never run.
  expect_error (scatterplot (iris [, -5], km$cluster, km$centers, ellipses = TRUE), NA)
  expect_error (scatterplot (iris [, -5], km$cluster, ellipses = TRUE), NA)
  # ... including when the cluster ids have a gap.
  expect_error (scatterplot (iris [, -5], ifelse (km$cluster == 2, 3, km$cluster),
                             ellipses = TRUE), NA)
})

# --- Regression: intern.intraclass() indexed the centers by cluster value ----------------

test_that ("intern.intraclass() does not depend on how the clusters are numbered", {
  data (iris)
  clus = c (rep (0, 50), rep (1, 50), rep (2, 50))
  expect_equal (intern.intraclass (clus, iris [, -5]),
                intern.intraclass (clus + 1, iris [, -5]))
  expect_error (intern.intraclass (c (rep (1, 50), rep (3, 100)), iris [, -5]), NA)
})

test_that ("intern.intraclass() agrees with kmeans' own within-cluster sum of squares", {
  data (iris)
  set.seed (1)
  km = stats::kmeans (iris [, -5], 3)
  expect_equal (intern.intraclass (km$cluster, iris [, -5]), km$tot.withinss)
})

# --- Regression: treeplot() ---------------------------------------------------------------

test_that ("treeplot() copes with k = 1, missing labels and a raw agnes object", {
  skip_if_not_installed ("cluster")
  data (iris)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  hca = HCA (iris [, -5], k = 3, method = "ward")
  expect_error (treeplot (hca), NA)
  expect_error (treeplot (hca, k = 1), NA)
  # A clustering with no labels: max (strwidth (NULL)) used to be -Inf, and par() then failed.
  expect_error (treeplot (HCA (as.matrix (unname (iris [, -5])), k = 3), labels = TRUE), NA)
  # The documentation says agnes results are accepted; rect.hclust() needs an hclust.
  expect_error (treeplot (cluster::agnes (iris [, -5], method = "ward"), k = 3), NA)
})

# ========================================================================================
# Regression tests for the "batch 5" fixes (performance).
# ========================================================================================

# --- Regression: comp = "pairwise" used to materialise one row per pair -------------------
# comparison.matrix() built an n (n - 1) / 2 by n matrix -- 12.8 MB for n = 150 and 29.8 GB
# for n = 2000 -- purely to count how the two partitions classify each pair. The counts now
# come from the contingency table, and the values must be identical.

test_that ("the pairwise counts agree with an explicit enumeration of the pairs", {
  set.seed (1)
  for (trial in 1:4)
  {
    n = 60
    k1 = sample (1:4, n, TRUE)
    k2 = if (trial == 3) k1 else if (trial == 4) ifelse (k1 <= 2, 1, 2) else sample (1:3, n, TRUE)
    pairs = utils::combn (n, 2)
    sep1 = k1 [pairs [1, ]] != k1 [pairs [2, ]]
    sep2 = k2 [pairs [1, ]] != k2 [pairs [2, ]]
    counts = fdm2id:::pairwise.counts (k1, k2)
    expect_equal (counts$both,    sum (!sep1 & !sep2), info = trial)
    expect_equal (counts$only1,   sum (!sep1 &  sep2), info = trial)
    expect_equal (counts$only2,   sum ( sep1 & !sep2), info = trial)
    expect_equal (counts$neither, sum ( sep1 &  sep2), info = trial)
    expect_equal (counts$total,   ncol (pairs), info = trial)
    # ... and so do the two indices built on them.
    expect_equal (fdm2id:::pairwise.rand (k1, k2), mean (sep1 == sep2), info = trial)
  }
})

test_that ("pairwise.kappa() reproduces irr::kappa2() on the explicit pair vectors", {
  skip_if_not_installed ("irr")
  set.seed (2)
  n = 60
  k1 = sample (1:4, n, TRUE)
  k2 = sample (1:3, n, TRUE)
  pairs = utils::combn (n, 2)
  sep1 = k1 [pairs [1, ]] != k1 [pairs [2, ]]
  sep2 = k2 [pairs [1, ]] != k2 [pairs [2, ]]
  expect_equal (fdm2id:::pairwise.kappa (k1, k2),
                irr::kappa2 (cbind (sep1, sep2), weight = "equal")$value)
})

test_that ("compare (comp = 'pairwise') is fast enough to be usable", {
  data (iris)
  set.seed (1)
  km = stats::kmeans (iris [, -5], 3)
  for (e in c ("accuracy", "jaccard", "kappa"))
    expect_true (is.finite (compare (km$cluster, iris [, 5], eval = e, comp = "pairwise")),
                 info = e)
  # 2000 observations used to need a 29.8 GB matrix; it must now be instantaneous.
  set.seed (3)
  big1 = sample (1:5, 2000, TRUE)
  big2 = sample (1:3, 2000, TRUE)
  expect_lt (system.time (fdm2id:::pairwise.rand (big1, big2)) [["elapsed"]], 1)
})

# =========================================================================================
# Sixth batch of the audit
# =========================================================================================

# --- Regression: DBSCAN() demanded both parameters with no way to guess them ---------------

test_that ("DBSCAN() has usable defaults and still honours explicit ones", {
  skip_if_not_installed ("fpc")
  data (iris)
  expect_equal (eval (formals (DBSCAN)$minpts), 5)
  expect_null (eval (formals (DBSCAN)$eps))
  auto = DBSCAN (iris [, -5])
  expect_true (max (auto$cluster) >= 1)
  # An explicit eps bypasses the estimation entirely.
  fixed = DBSCAN (iris [, -5], minpts = 5, eps = 0.5)
  expect_equal (fixed$cluster, DBSCAN (iris [, -5], minpts = 5, eps = 0.5)$cluster)
})

test_that ("dbscan.eps() returns a positive distance inside the data's range", {
  data (iris)
  eps = fdm2id:::dbscan.eps (iris [, -5], 5)
  expect_true (eps > 0)
  expect_true (eps < max (stats::dist (iris [, -5])))
})

# --- Regression: KMEANS() had no seed, although stats::kmeans (nstart = 10) is random -------

test_that ("KMEANS() is reproducible when seeded", {
  data (iris)
  expect_true ("seed" %in% names (formals (KMEANS)))
  a = KMEANS (iris [, -5], k = 3, seed = 0)
  b = KMEANS (iris [, -5], k = 3, seed = 0)
  expect_equal (a$cluster, b$cluster)
})

# --- Regression: 2:max ran backwards for max < 2 -------------------------------------------

test_that ("kmeans.getk() refuses a search range it cannot search", {
  data (iris)
  expect_error (fdm2id:::kmeans.getk (iris [, -5], max = 1), "at least 2")
  expect_error (fdm2id:::kmeans.getk (iris [, -5], max = 3, graph = FALSE, seed = 0), NA)
})

# =========================================================================================
# Eighth batch: the three pairwise comparison indices
# =========================================================================================

# --- Regression: compare.jaccard (comp = "pairwise") returned the Rand index ---------------

test_that ("the three pairwise criteria are three different indices", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, seed = 0, graph = FALSE)
  rand = compare.accuracy (km$cluster, iris [, 5], comp = "pairwise")
  jacc = compare.jaccard (km$cluster, iris [, 5], comp = "pairwise")
  ari = compare.kappa (km$cluster, iris [, 5], comp = "pairwise")
  expect_false (isTRUE (all.equal (rand, jacc)))
  expect_false (isTRUE (all.equal (rand, ari)))
  expect_true (jacc < rand)
})

test_that ("the pairwise indices match their definitions on hand-computed counts", {
  # 4 observations, 6 pairs. k1 groups {1,2} and {3,4}; k2 groups {1,2,3} and {4}.
  k1 = c (1, 1, 2, 2)
  k2 = c (1, 1, 1, 2)
  p = fdm2id:::pairwise.counts (k1, k2)
  expect_equal (p$total, 6)
  expect_equal (p$both, 1)      # la paire (1,2)
  expect_equal (p$only1, 1)     # (3,4)
  expect_equal (p$only2, 2)     # (1,3) et (2,3)
  expect_equal (p$neither, 2)   # (1,4) et (2,4)
  expect_equal (fdm2id:::pairwise.rand (k1, k2), (1 + 2) / 6)
  expect_equal (fdm2id:::pairwise.jaccard (k1, k2), 1 / (1 + 1 + 2))
})

test_that ("the pairwise indices behave at the extremes", {
  data (iris)
  gt = iris [, 5]
  n = nrow (iris)
  # Identical partitions: everything agrees.
  expect_equal (compare.jaccard (as.numeric (gt), gt, comp = "pairwise"), 1)
  expect_equal (compare.accuracy (as.numeric (gt), gt, comp = "pairwise"), 1)
  # One observation per cluster: not a single pair is grouped together, so the Jaccard index
  # is 0 -- while the Rand index stays high because every pair is (trivially) kept apart.
  expect_equal (compare.jaccard (1:n, gt, comp = "pairwise"), 0)
  expect_gt (compare.accuracy (1:n, gt, comp = "pairwise"), 0.6)
  # Everything in one cluster: no pair is separated, so Rand and Jaccard coincide.
  expect_equal (compare.jaccard (rep (1, n), gt, comp = "pairwise"),
                compare.accuracy (rep (1, n), gt, comp = "pairwise"))
})

test_that ("compare.kappa (comp = 'pairwise') is Cohen's kappa on the pair table", {
  skip_if_not_installed ("irr")
  data (iris)
  gt = as.numeric (iris [, 5])
  pairs = utils::combn (150, 2)
  cases = list (KMEANS (iris [, -5], k = 3, seed = 0, graph = FALSE)$cluster,
                KMEANS (iris [, -5], k = 10, seed = 0, graph = FALSE)$cluster,
                sample (1:4, 150, TRUE))
  for (k in cases)
  {
    together = cbind (as.numeric (k [pairs [1, ]] == k [pairs [2, ]]),
                      as.numeric (gt [pairs [1, ]] == gt [pairs [2, ]]))
    expect_equal (compare.kappa (k, iris [, 5], comp = "pairwise"),
                  irr::kappa2 (together)$value)
  }
})

test_that ("Cohen's kappa on the pair table coincides with the adjusted Rand index", {
  # Warrens (2008), Journal of Classification 25(2):177-183. This is a property of the
  # statistic, not the definition compare.kappa() implements -- see the test above.
  skip_if_not_installed ("mclust")
  data (iris)
  gt = as.numeric (iris [, 5])
  cases = list (KMEANS (iris [, -5], k = 3, seed = 0, graph = FALSE)$cluster,
                KMEANS (iris [, -5], k = 10, seed = 0, graph = FALSE)$cluster,
                gt, rep (1, 150), 1:150)
  for (k in cases)
    expect_equal (compare.kappa (k, iris [, 5], comp = "pairwise"),
                  mclust::adjustedRandIndex (k, gt))
})

test_that ("the other comparison modes were already distinct and stay so", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, seed = 0, graph = FALSE)
  for (m in c ("max", "cluster"))
    expect_false (isTRUE (all.equal (compare.accuracy (km$cluster, iris [, 5], comp = m),
                                     compare.jaccard (km$cluster, iris [, 5], comp = m))))
})

# =========================================================================================
# Ninth batch: choosing the number of clusters
# =========================================================================================

test_that ("the four criteria all return a usable number of clusters on iris", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  for (crit in c ("pseudo-F", "silhouette", "elbow"))
  {
    k = kmeans.getk (iris [, -5], criterion = crit, seed = 0)
    expect_true (k >= 2 && k <= 9, info = crit)
  }
  skip_if_not_installed ("cluster")
  k = kmeans.getk (iris [, -5], criterion = "gap", B = 10, seed = 0)
  expect_true (k >= 1 && k <= 9)
})

test_that ("KMEANS() honours the criterion it is given", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  for (crit in c ("pseudo-F", "silhouette", "elbow"))
  {
    km = KMEANS (iris [, -5], k = 9, criterion = crit, seed = 0)
    expect_equal (length (km$size), kmeans.getk (iris [, -5], max = 9, criterion = crit,
                                                 graph = FALSE, seed = 0), info = crit)
  }
  # criterion = "none" uses k as it is.
  expect_equal (length (KMEANS (iris [, -5], k = 4, seed = 0)$size), 4)
  expect_error (KMEANS (iris [, -5], criterion = "banana"), "should be one of")
})

test_that ("the gap statistic can answer 'no cluster structure at all'", {
  skip_if_not_installed ("cluster")
  set.seed (0)
  noise = matrix (stats::rnorm (400), ncol = 2)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # This is what the gap statistic is for, and what none of the other three can say: they are
  # all maximised over k >= 2 and must return a number of clusters whatever the data.
  expect_equal (kmeans.getk (noise, max = 6, criterion = "gap", B = 30, seed = 0), 1)
  expect_gte (kmeans.getk (noise, max = 6, criterion = "silhouette", seed = 0), 2)
})

test_that ("silhouette and pseudo-F agree with a direct computation", {
  skip_if_not_installed ("cluster")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  dis = stats::dist (iris [, -5])
  direct = sapply (2:6, function (k)
  {
    set.seed (0)
    mean (cluster::silhouette (stats::kmeans (iris [, -5], k, nstart = 10)$cluster,
                               dis) [, "sil_width"])
  })
  expect_equal (kmeans.getk (iris [, -5], max = 6, criterion = "silhouette", graph = FALSE,
                             seed = 0),
                (2:6) [which.max (direct)])
})

test_that ("knee.index() finds the bend of a curve, and copes with degenerate ones", {
  # Two straight segments meeting at the 4th point.
  y = c (10, 7, 4, 1, .9, .8, .7)
  expect_equal (fdm2id:::knee.index (y), 4)
  expect_equal (fdm2id:::knee.index (rep (1, 10)), 5)   # flat: no bend to find
  expect_equal (fdm2id:::knee.index (c (1, 2)), 1)      # too short
})

# --- The print methods of the clustering results --------------------------------------------

test_that ("the clustering print methods say what was found", {
  data (iris)
  shown = function (x) paste (utils::capture.output (print (x)), collapse = "\n")
  out = shown (DBSCAN (iris [, -5], minpts = 5, eps = 0.65))
  expect_match (out, "DBSCAN clustering")
  expect_match (out, "unclustered")
  expect_match (shown (EM (iris [, -5], 3)), "EM clustering")
  expect_match (shown (EM (iris [, -5], 3)), "log-likelihood")
  skip_if_not_installed ("kohonen")
  expect_match (shown (SOM (iris [, -5], 4, 4)), "4 x 4")
})

# =========================================================================================
# Tenth batch: PAM
# =========================================================================================

test_that ("PAM() partitions around actual observations of the dataset", {
  skip_if_not_installed ("cluster")
  data (iris)
  model = PAM (iris [, -5], 3)
  expect_equal (length (unique (model$cluster)), 3)
  expect_equal (nrow (model$medoids), 3)
  # A medoid is an observation of the dataset, not an average of them.
  for (i in 1:3)
    expect_true (any (apply (iris [, -5], 1, function (r) isTRUE (all.equal (unname (r),
                                                                            unname (model$medoids [i, ]))))))
  expect_equal (unname (model$centers), unname (model$medoids))
})

test_that ("PAM() can choose k by silhouette, and refuses an impossible range", {
  skip_if_not_installed ("cluster")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  k = length (unique (PAM (iris [, -5], k = 6, criterion = "silhouette", seed = 0)$cluster))
  expect_true (k >= 2 && k <= 6)
  expect_error (PAM (iris [, -5], k = 1, criterion = "silhouette"), "at least 2")
})

test_that ("the new clustering results can be plotted like the others", {
  skip_if_not_installed ("cluster")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (plotclus (PAM (iris [, -5], 3), iris [, -5]), NA)
})

# =========================================================================================
# Eleventh batch
# =========================================================================================

test_that ("scatterplot() accepts the factors its documentation promises", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # The body treated 'clusters' as numeric throughout, so a factor stopped it on
  # "'min' not meaningful for factors" and a character vector on "non-numeric argument to
  # binary operator" -- while the documentation announced both.
  expect_error (scatterplot (iris [, -5], iris [, 5]), NA)
  expect_error (scatterplot (iris [, -5], as.character (iris [, 5])), NA)
  expect_error (scatterplot (iris [, -5], as.numeric (iris [, 5])), NA)
  expect_error (scatterplot (iris [, -5], iris [, 5], ellipses = TRUE), NA)
  # A density clustering marks its noise as cluster 0, which must keep working.
  db = suppressMessages (DBSCAN (iris [, -5]))
  expect_error (scatterplot (iris [, -5], db$cluster), NA)
})

test_that ("kdistances() agrees with sorting every column, and is not O(n log n) per column", {
  data (iris)
  reference = apply (as.matrix (stats::dist (iris [, -5])), 2, sort) [6, ]
  expect_equal (unname (fdm2id:::kdistances (iris [, -5], 5)), unname (reference))
  expect_equal (unname (fdm2id:::kdistances (iris [, -5], 1)),
                unname (apply (as.matrix (stats::dist (iris [, -5])), 2, sort) [2, ]))
})

# =========================================================================================
# Third audit, D5
# =========================================================================================

test_that ("boxclus accepts the class labels its documentation promises", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  km = KMEANS (iris [, -5], k = 3, seed = 0)
  expect_error (boxclus (iris [, -5], km$cluster), NA)
  # A factor used to fail on "'min' not meaningful for factors", where scatterplot() -- the
  # other plot of the same clusters -- accepted it.
  expect_error (boxclus (iris [, -5], iris [, 5]), NA)
  expect_error (boxclus (iris [, -5], as.character (iris [, 5])), NA)
  db = suppressMessages (DBSCAN (iris [, -5], minpts = 5, eps = .4))
  expect_error (boxclus (iris [, -5], db$cluster), NA)
})

test_that ("the clustering legend names the classes, the noise and their colours", {
  data (iris)
  # Colours follow the 1 + code convention, so the legend matches the boxes and the points.
  # boxclus() used to colour its boxes 2:(k + 1) while reading the legend off 1 + clusters,
  # which shifted every entry as soon as there was noise.
  noise = c (rep (0, 5), rep (1, 5), rep (2, 5))
  expect_equal (fdm2id:::cluster.legend (noise),
                list (labels = c ("Noise", "Cluster 1", "Cluster 2"), col = c (1, 2, 3)))
  expect_equal (fdm2id:::cluster.legend (c (1, 1, 2, 2)),
                list (labels = c ("Cluster 1", "Cluster 2"), col = c (2, 3)))
  coded = fdm2id:::cluster.codes (iris [, 5])
  expect_equal (coded$names, levels (iris [, 5]))
  expect_true (is.numeric (coded$codes))
  expect_equal (fdm2id:::cluster.legend (coded$codes, coded$names)$labels, levels (iris [, 5]))
  # A numeric vector is left exactly as it is.
  expect_equal (fdm2id:::cluster.codes (noise), list (codes = noise, names = NULL))
})

# =========================================================================================
# Fourth audit, B7 / B11 / D10
# =========================================================================================

test_that ("Dunn's index refuses a partition it cannot be computed on", {
  data (iris)
  # A single cluster has nothing to be separated from: the inner min() ran over an empty set,
  # warned about it, and returned Inf.
  expect_error (intern.dunn (rep (1, nrow (iris)), iris [, -5]), "at least two")
  expect_error (intern (rep (1, nrow (iris)), iris [, -5], eval = "dunn"), "at least two")
  # Every cluster a singleton: all diameters are zero, so the index is 0 / 0.
  x = data.frame (a = 1:5, b = 1:5)
  expect_equal (intern.dunn (1:5, x), Inf)
  expect_true (all (is.infinite (intern.dunn (1:5, x, type = "cluster"))))
})

test_that ("the per-cluster Dunn values decompose the global one", {
  data (iris)
  # Both divide by the largest diameter of the partition, so the global index is the minimum
  # of the per-cluster ones. Each used to divide by its own cluster's diameter, which made them
  # two different statistics wearing the same name.
  for (k in 2:4)
  {
    clus = KMEANS (iris [, -5], k = k, seed = 0)$cluster
    parts = intern.dunn (clus, iris [, -5], type = "cluster")
    expect_length (parts, k)
    expect_equal (min (parts), intern.dunn (clus, iris [, -5]))
  }
  # A singleton no longer sends its own value to infinity.
  clus = c (rep (1, 74), rep (2, 75), 3)
  parts = intern.dunn (clus, iris [, -5], type = "cluster")
  expect_true (all (is.finite (parts)))
  expect_equal (min (parts), intern.dunn (clus, iris [, -5]))
})

test_that ("the clustering criteria name a misspelt 'type' or 'comp'", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, seed = 0)
  # These used to be read as 'type [1] == "global"' / 'comp [1] == "pairwise"', so anything
  # else silently answered a different question.
  for (f in list (intern.dunn, intern.interclass, intern.intraclass))
    expect_error (f (km$cluster, iris [, -5], type = "clsuter"), "should be one of")
  expect_error (intern (km$cluster, iris [, -5], type = "clsuter"), "should be one of")
  for (f in list (compare.accuracy, compare.jaccard, compare.kappa))
    expect_error (f (km$cluster, iris [, 5], comp = "pairwize"), "should be one of")
  expect_error (compare (km$cluster, iris [, 5], comp = "pairwize"), "should be one of")
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (plotclus (km, iris [, -5], type = "scatterplot"), "should be one of")
  # A valid type that this clustering cannot draw still says so, by name.
  expect_error (plotclus (km, iris [, -5], type = "tree"), "hierarchical")
})

test_that ("spectral clustering diagonalises a symmetric matrix and says when it cannot", {
  data (iris)
  res = SPECTRAL (iris [, -5], 3, seed = 0)
  expect_length (res$cluster, nrow (iris))
  expect_false (any (is.complex (res$proj)))
  expect_equal (compare (res$cluster, iris [, 5]) > 0.7, c (accuracy = TRUE))
  # A sigma far too small leaves every affinity at zero, and the normalisation divided by it.
  expect_error (SPECTRAL (iris [, -5], 3, sigma = 1e-8), "too small")
})

# =========================================================================================
# Sixth audit: cost of the hand-written algorithms, at identical results
# =========================================================================================

test_that ("the k leading eigenvectors give the clustering the full decomposition gives", {
  data (iris)
  # SPECTRAL asks RSpectra for the k eigenvectors it needs instead of all n of them (143 s
  # down to 3 s on 5000 observations). Eigenvectors are defined up to a sign, and the
  # projection is row-normalised before k-means, so the partition must be the same either way.
  a = exp (-flexclust::dist2 (iris [, -5], iris [, -5])^2 / 2)
  diag (a) = 0
  p = 1 / sqrt (rowSums (a))
  l = a * outer (p, p)
  full = eigen (l, symmetric = TRUE)$vectors [, 1:3]
  quick = fdm2id:::spectral.eigenvectors (l, 3)
  partition = function (x)
  {
    proj = sweep (x, 1, sqrt (rowSums (x * x)), "/")
    set.seed (0)
    stats::kmeans (proj, centers = 3, nstart = 100)$cluster
  }
  expect_equal (mclust::adjustedRandIndex (partition (full), partition (quick)), 1)
  # The fallback is exercised whenever RSpectra is not installed, and must still work.
  expect_equal (dim (quick), c (nrow (iris), 3))
})

test_that ("Dunn's index is computed block by block and still matches fpc", {
  skip_if_not_installed ("fpc")
  data (iris)
  # The n x n distance matrix is no longer materialised in full (300 Mo down to 66 on 5000
  # observations); the numbers must not move.
  for (k in 2:4)
  {
    clus = KMEANS (iris [, -5], k = k, seed = 0)$cluster
    expect_equal (intern.dunn (clus, iris [, -5]),
                  fpc::cluster.stats (stats::dist (iris [, -5]), clus)$dunn)
  }
})

# =========================================================================================
# Seventh audit: HCA's engine
# =========================================================================================

test_that ("the two HCA engines build the same hierarchy", {
  skip_if_not_installed ("cluster")
  data (iris)
  # stats::hclust() is the default because it is far faster (0.6 s against 30 s on 3000
  # observations); it has to give the same answer, and "ward" maps to hclust's "ward.D2".
  for (method in c ("ward", "single", "average", "complete"))
  {
    quick = HCA (iris [, -5], k = 3, method = method)
    slow = HCA (iris [, -5], k = 3, method = method, engine = "agnes")
    expect_equal (mclust::adjustedRandIndex (quick$cluster, slow$cluster), 1, info = method)
    if (method != "complete") # complete-linkage ties are broken differently, same cut
      expect_equal (sort (quick$height), sort (slow$height), info = method)
  }
  # And the number of clusters the largest drop in height picks out.
  expect_equal (HCA (iris [, -5])$k, HCA (iris [, -5], engine = "agnes")$k)
})

test_that ("everything built on an HCA still works whichever engine produced it", {
  skip_if_not_installed ("cluster")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  for (engine in c ("hclust", "agnes"))
  {
    h = HCA (iris [, -5], k = 3, engine = engine)
    expect_s3_class (h, "hca")
    expect_true (inherits (h, "hclust"))
    expect_error (treeplot (h), NA, info = engine)
    expect_error (treeplot (h, labels = TRUE, k = 3, split = TRUE), NA, info = engine)
    expect_error (treeplot (h, horiz = TRUE), NA, info = engine)
    expect_error (plotclus (h, iris [, -5], type = "tree"), NA, info = engine)
    expect_error (plotclus (h, iris [, -5], type = "height"), NA, info = engine)
    expect_error (plotclus (h, iris [, -5], type = "scatter"), NA, info = engine)
    expect_length (predict (h, iris [, -5]), nrow (iris))
    expect_length (stats::cutree (h, 4), nrow (iris))
  }
})

test_that ("HCA names the linkages it accepts", {
  data (iris)
  expect_error (HCA (iris [, -5], k = 3, method = "flexible"), "unknown linkage")
  expect_error (HCA (iris [, -5], k = 3, method = "ward.D2"), NA)
  expect_error (HCA (iris [, -5], k = 3, engine = "hclast"), "should be one of")
})

test_that ("the Jaccard index is offered wherever a partition is compared", {
  data (iris)
  # accuracy compares two *binary indicators*, so the observations in neither the cluster nor
  # the class dominate the count and a random partition still scores about 0.6. The Jaccard
  # index does not count them, and has to stay available next to it.
  km = KMEANS (iris [, -5], k = 3, seed = 0)
  expect_true ("jaccard" %in% names (fdm2id:::compare.functions ()))
  expect_true ("jaccard" %in% names (fdm2id:::eval.functions ()))
  expect_equal (formals (stability)$eval, "jaccard")
  for (comp in c ("max", "cluster", "pairwise"))
    expect_true (all (is.finite (compare (km$cluster, iris [, 5], eval = "jaccard",
                                          comp = comp))), info = comp)
  expect_named (compare (km$cluster, iris [, 5], eval = c ("accuracy", "jaccard", "kappa")),
                c ("accuracy", "jaccard", "kappa"))
})

# =========================================================================================
# Tenth audit: stability when the number of clusters varies
# =========================================================================================

test_that ("stability copes with a method that finds its own number of clusters", {
  skip_if_not_installed ("fpc")
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  d = data.twomoons (seed = 0)
  # A bootstrap sample need not contain every cluster, and DBSCAN does not find the same
  # number twice: apply() then returned a list, is.vector() is TRUE for a list, and the whole
  # thing ended on mean (list) = NA with "argument is not numeric or logical".
  res = stability (DBSCAN, d [, -3], minpts = 4, eps = .2, seed = 0, nsampling = 5)
  expect_true (is.matrix (res))
  expect_true (all (is.finite (res)))
  expect_equal (colnames (res), "jaccard")
  global = stability (DBSCAN, d [, -3], minpts = 4, eps = .2, seed = 0, nsampling = 5,
                      type = "global")
  expect_named (global, "jaccard")
  expect_true (is.finite (global))
  # Several criteria at once, still per cluster.
  two = stability (DBSCAN, d [, -3], minpts = 4, eps = .2, seed = 0, nsampling = 5,
                   eval = c ("jaccard", "accuracy"))
  expect_equal (colnames (two), c ("jaccard", "accuracy"))
  expect_true (all (is.finite (two)))
})

test_that ("compare (comp = 'cluster') keeps the cluster names", {
  data (iris)
  km = KMEANS (iris [, -5], k = 3, seed = 0)
  # names (res) = eval overwrote the per-cluster names and padded with NA: two criteria on
  # three clusters came back as six values named "jaccard", "accuracy", NA, NA, NA, NA.
  one = compare (km$cluster, iris [, 5], eval = "jaccard", comp = "cluster")
  expect_named (one, paste ("Cluster", 1:3))
  two = compare (km$cluster, iris [, 5], eval = c ("jaccard", "accuracy"), comp = "cluster")
  expect_true (is.matrix (two))
  expect_equal (rownames (two), c ("jaccard", "accuracy"))
  expect_equal (colnames (two), paste ("Cluster", 1:3))
  # The other two modes keep answering one number per criterion.
  expect_named (compare (km$cluster, iris [, 5], eval = c ("accuracy", "jaccard")),
                c ("accuracy", "jaccard"))
  expect_named (compare (km$cluster, iris [, 5], eval = c ("accuracy", "jaccard"),
                         comp = "pairwise"), c ("accuracy", "jaccard"))
})

# --- Regression: a character ground truth was coerced with as.numeric() ------------------
# read.table() stopped converting strings to factors in R 4.0, so gt arrives as a character
# vector where it used to arrive as a factor. as.numeric() then gave NA throughout and the
# three indices came back as a plausible-looking 1, with only a coercion warning.

test_that ("compare.* give the same value for a factor and a character ground truth", {
  data (iris)
  clus = as.numeric (iris [, 5])
  clus [c (1, 60, 120)] = c (2, 3, 1) # not a perfect partition, so 1 is not the right answer
  for (comp in c ("max", "pairwise"))
    for (fun in list (compare.accuracy, compare.jaccard, compare.kappa))
    {
      ref = expect_silent (fun (clus, iris [, 5], comp = comp))
      expect_equal (fun (clus, as.character (iris [, 5]), comp = comp), ref)
      expect_lt (unname (ref), 1)
    }
})
