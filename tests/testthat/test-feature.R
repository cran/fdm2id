# Tests for feature.R

# --- Regression: selectfeatures()'s fs.*/fseval.* dispatch tables (get(paste(...)) fix) --

test_that ("selectfeatures() ranking algorithm works and respects uninb", {
  data (iris)
  res = selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking", uninb = 2)
  expect_equal (length (res$selection), 2)
  expect_s3_class (res, "selection")
})

test_that ("selectfeatures() forward algorithm works with a multivariate criterion", {
  data (iris)
  res = selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "fstat")
  expect_true (length (res$selection) >= 1)
  expect_true (length (res$selection) <= ncol (iris) - 1)
})

test_that ("selectfeatures() rejects an unknown algorithm/criterion with a clear error", {
  data (iris)
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "not_a_real_algorithm"))
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
                                unieval = "not_a_real_criterion"))
})

test_that ("FEATURESELECTION() smoke test: selection + downstream model + predict", {
  skip_if_not_installed ("MASS")
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  model = FEATURESELECTION (d$train.x, d$train.y, uninb = 2, mainmethod = LDA)
  pred = predict (model, d$test.x)
  expect_equal (length (pred), nrow (d$test.x))
})

# --- Regression: selectfeatures() used to return a malformed "selection" object ----------
# fs.ranking() message()d and returned NULL when it had no criterion to work with;
# selectfeatures() then wrote res$univariate <- ... on that NULL, which *creates* a list
# instead of failing, and handed back a "selection" object holding nothing but a name.

test_that ("selectfeatures() says what is missing instead of returning an empty object", {
  data (iris)
  # 'ranking' with no criterion at all used to hand back an empty "selection" object; it then
  # raised an error, and since the uniformity review it simply works -- 'multieval' now has a
  # default there too, which chooses how many of the ranked features to keep.
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking"), NA)
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
                                unieval = "fisher", uninb = NULL, unithreshold = NULL,
                                multieval = NULL),
                "uninb")
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
                                multieval = "wrapper"),
                "wrapmethod")
})

test_that ("the documented wrapper example of selectfeatures() runs", {
  skip_if_not_installed ("MASS")
  skip_if_not_installed ("e1071")
  data (iris)
  res = selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
                        multieval = "wrapper", wrapmethod = LDA)
  expect_s3_class (res, "selection")
  expect_gt (length (res$selection), 0)
})

test_that ("selectfeatures (keep = TRUE) keeps a data.frame even for a single feature", {
  data (iris)
  res = selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking", uninb = 1, keep = TRUE)
  expect_equal (ncol (res$dataset), 1)
  expect_equal (nrow (res$dataset), nrow (iris))
})

# --- Regression: mutual information was estimated by kernel density -----------------------
# proba() rescanned the whole sample for every cell of a Sturges grid, so mRMR needed 4
# seconds on iris (150 x 4), and the estimate was not a mutual information: it fed
# stats::cor() into a formula expecting a covariance matrix. It was first replaced by a
# discretised estimator, then (eighth batch, see below) by the nearest-neighbour estimators
# that are the reference for this problem. Values differ from previous versions.

test_that ("mutual information has the properties it should", {
  data (iris)
  set.seed (1)
  x = stats::rnorm (500)
  y = stats::rnorm (500)
  # Non-negative, symmetric, near zero for independent variables, and much larger for a
  # variable against itself.
  expect_gt (fdm2id:::mutualinformation (x, x), fdm2id:::mutualinformation (x, y))
  expect_lt (fdm2id:::mutualinformation (x, y), 0.5)
  expect_equal (fdm2id:::mutualinformation (iris [, 1], iris [, 3]),
                fdm2id:::mutualinformation (iris [, 3], iris [, 1]))
  # Petal length says more about the species than sepal width does.
  expect_gt (fdm2id:::mutualinformation (iris [, 3], iris [, 5]),
             fdm2id:::mutualinformation (iris [, 2], iris [, 5]))
  # A constant variable carries no information.
  expect_equal (fdm2id:::mutualinformation (rep (1, 100), stats::rnorm (100)), 0)
  # Affine invariance is checked on a genuine pair of variables (see the eighth-batch tests).
  # It is deliberately *not* checked on the degenerate mutualinformation (x, x): the mutual
  # information of a continuous variable with itself is infinite, so what the estimator
  # returns there is an artefact of the sample, not a quantity worth pinning down.
  expect_equal (fdm2id:::mutualinformation (x, y), fdm2id:::mutualinformation (2 * x + 1, y))
})

test_that ("selectfeatures() with mrmr is usable", {
  data (iris)
  expect_lt (system.time (res <- selectfeatures (iris [, -5], iris [, 5],
                                                 algorithm = "forward",
                                                 multieval = "mrmr")) [["elapsed"]], 2)
  expect_gt (length (res$selection), 0)
})

# =========================================================================================
# Sixth batch of the audit
# =========================================================================================

# --- Regression: FEATURESELECTION() and selectfeatures() disagreed on their defaults --------

test_that ("both feature-selection entry points share one list of criteria", {
  defaults = function (f, argument)
  {
    e = new.env ()
    assign ("algorithm", "ranking", envir = e)
    return (eval (formals (f) [[argument]], envir = e))
  }
  expect_equal (defaults (FEATURESELECTION, "unieval"), defaults (selectfeatures, "unieval"))
  expect_equal (defaults (FEATURESELECTION, "multieval"),
                defaults (selectfeatures, "multieval"))
  expect_true ("mrmr" %in% fdm2id:::fseval.multivariate ())
  # Every advertised criterion is actually implemented.
  known = names (fdm2id:::fseval.functions ())
  expect_true (all (fdm2id:::fseval.univariate () %in% known))
  expect_true (all (setdiff (fdm2id:::fseval.multivariate (), "wrapper") %in% known))
})

# =========================================================================================
# Eighth batch: the mutual information behind mRMR
# =========================================================================================

# --- The properties a mutual information must have. The 0.9.10 estimator failed the first
# --- three of them, and the histogram estimator that briefly replaced it was heavily biased.

test_that ("mutualinformation() is never negative", {
  mi = fdm2id:::mutualinformation
  set.seed (0)
  expect_true (all (replicate (20, mi (stats::rnorm (150), stats::rnorm (150))) >= 0))
  expect_true (all (replicate (20, mi (stats::rnorm (150),
                                       factor (sample (1:3, 150, TRUE)))) >= 0))
  expect_gte (mi (factor (sample (1:3, 150, TRUE)), factor (sample (1:2, 150, TRUE))), 0)
})

test_that ("mutualinformation() does not depend on the units of the variables", {
  mi = fdm2id:::mutualinformation
  set.seed (0)
  x = stats::rnorm (300)
  y = .6 * x + sqrt (1 - .36) * stats::rnorm (300)
  ref = mi (x, y)
  expect_equal (mi (1000 * x, y), ref)
  expect_equal (mi (x, 1000 * y), ref)
  expect_equal (mi (x + 50, y), ref)
  k = factor (sample (c ("a", "b"), 300, TRUE))
  expect_equal (mi (1000 * x, k), mi (x, k))
})

test_that ("mutualinformation() is symmetric", {
  mi = fdm2id:::mutualinformation
  set.seed (0)
  x = stats::rnorm (200)
  y = .5 * x + stats::rnorm (200)
  k = factor (sample (1:3, 200, TRUE))
  expect_equal (mi (x, y), mi (y, x))
  expect_equal (mi (x, k), mi (k, x))
})

test_that ("mutualinformation() recovers the known value of a bivariate normal", {
  mi = fdm2id:::mutualinformation
  set.seed (1)
  n = 2000
  for (rho in c (0, .6, .9))
  {
    x = stats::rnorm (n)
    y = rho * x + sqrt (1 - rho^2) * stats::rnorm (n)
    truth = -0.5 * log2 (1 - rho^2)
    # The estimator is consistent, not exact: 0.05 bit is well inside its sampling error and
    # far below the 0.3 bit bias a binned estimator shows on this problem.
    expect_lt (abs (mi (x, y) - truth), 0.05)
  }
})

test_that ("mutualinformation() is deterministic and leaves the user's RNG untouched", {
  mi = fdm2id:::mutualinformation
  set.seed (0)
  x = round (stats::rnorm (200), 1)   # plenty of ties, so the tie-breaking runs
  y = round (stats::rnorm (200), 1)
  expect_identical (mi (x, y), mi (x, y))
  set.seed (123)
  before = stats::runif (1)
  set.seed (123)
  invisible (mi (x, y))
  expect_identical (stats::runif (1), before)
})

test_that ("mutualinformation() checks that its two arguments match", {
  expect_error (fdm2id:::mutualinformation (1:10, 1:9), "same length")
})

# --- Regression: the estimator returned negative values on the package's own data ----------

test_that ("no variable of wine or spine gets a negative mutual information", {
  data (wine)
  data (spine)
  expect_true (all (sapply (wine [, -1],
                            function (v) fdm2id:::mutualinformation (v, wine [, 1])) >= 0))
  expect_true (all (sapply (spine [, 1:6],
                            function (v) fdm2id:::mutualinformation (v, spine$Classif3)) >= 0))
})

# --- The cache must not change any result --------------------------------------------------

test_that ("fseval.mrmr() returns the same score with and without its cache", {
  data (iris)
  x = iris [, -5]
  y = iris [, 5]
  direct = fdm2id:::fseval.mrmr (x, y, vtype = "multivariate")
  cache = new.env (parent = emptyenv ())
  first = fdm2id:::fseval.mrmr (x, y, vtype = "multivariate", micache = cache)
  second = fdm2id:::fseval.mrmr (x, y, vtype = "multivariate", micache = cache)
  expect_equal (first, direct)
  expect_equal (second, direct)
  # 4 relevance terms + 6 pairs, each stored once.
  expect_equal (length (ls (cache)), 10)
})

test_that ("fseval.mrmr() keeps the types of a mixed data.frame", {
  data (iris)
  mixed = data.frame (num = iris [, 1], fac = iris [, 5])
  # as.matrix() on such a data.frame yields a character matrix, which would make the numeric
  # column be estimated as if it were categorical.
  expect_false (is.na (fdm2id:::fseval.mrmr (mixed, iris [, 5], vtype = "multivariate")))
})

test_that ("selectfeatures (multieval = 'mrmr') runs on the usual datasets", {
  data (iris)
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "forward",
                                multieval = "mrmr"), NA)
  expect_error (selectfeatures (iris [, -5], iris [, 5], algorithm = "backward",
                                multieval = "mrmr"), NA)
  # A single feature has no redundancy term, only relevance.
  expect_gte (fdm2id:::fseval.mrmr (iris [, 1], iris [, 5], vtype = "multivariate"), 0)
})

# =========================================================================================
# Eleventh batch: Relief rewritten
# =========================================================================================

# A slow but indisputable Relief-F, straight from the definition (Kononenko, 1994), used to
# check the vectorised implementation rather than checking it against its own past output.
relief.reference <-
  function (train, labels, samples, k = 10)
  {
    train = as.matrix (train)
    labels = factor (labels)
    n = nrow (train)
    lev = levels (labels)
    prior = as.vector (table (labels)) / n
    names (prior) = lev
    diffv = apply (train, 2, function (v) diff (range (v)))
    diffv [diffv == 0] = 1
    w = numeric (ncol (train))
    for (s in samples)
      for (l in lev)
      {
        cls = which (labels == l)
        d = sqrt (colSums ((t (train [cls, , drop = FALSE]) - train [s, ])^2))
        o = cls [order (d)]
        if (l == labels [s])
          o = o [-1]
        kl = min (k, length (o))
        if (kl < 1)
          next
        neighbours = o [1:kl]
        diff = sweep (abs (sweep (train [neighbours, , drop = FALSE], 2, train [s, ], "-")),
                      2, diffv, "/")
        coefficient = if (l == labels [s]) -1
                      else prior [l] / (1 - prior [as.character (labels [s])])
        w = w + coefficient * colSums (diff) / (length (samples) * kl)
      }
    return (w)
  }

test_that ("fseval.relief() computes what Relief-F is defined to be", {
  data (iris)
  set.seed (1)
  a = fdm2id:::fseval.relief (iris [, -5], iris [, 5])
  set.seed (1)
  b = relief.reference (iris [, -5], iris [, 5], sample (nrow (iris), nrow (iris)))
  expect_equal (unname (a), unname (b))
})

test_that ("fseval.relief() no longer returns NA when a class is smaller than k", {
  data (iris)
  # 50 / 6 / 50: the second class has fewer than the k + 1 = 11 neighbours the old
  # 'order (v) [1:(k + 1)]' asked for, so it produced NA indices and every weight came back
  # as NA -- silently, and with it every feature selection built on Relief or CFS.
  d = rbind (iris [1:50, ], iris [51:56, ], iris [101:150, ])
  d$Species = droplevels (d$Species)
  set.seed (1)
  a = fdm2id:::fseval.relief (d [, -5], d [, 5])
  expect_false (any (is.na (a)))
  set.seed (1)
  expect_equal (unname (a),
                unname (relief.reference (d [, -5], d [, 5], sample (nrow (d), nrow (d)))))
  expect_false (any (is.na (fdm2id:::fseval.cfs (d [, -5], d [, 5]))))
})

test_that ("fseval.relief() normalises each variable by its own range", {
  # 'abs (...) / diffv' on a (k x p) matrix recycles down the columns, so each value was
  # divided by some other variable's range. Making one variable a thousand times larger than
  # the others must not change the ranking, since Relief divides by the range.
  data (iris)
  set.seed (1)
  a = fdm2id:::fseval.relief (iris [, -5], iris [, 5])
  scaled = iris [, -5]
  scaled [, 1] = scaled [, 1] * 1000
  set.seed (1)
  b = fdm2id:::fseval.relief (scaled, iris [, 5])
  expect_equal (order (a), order (b))
})

test_that ("knn.indices() returns the k nearest, dropping the point itself when asked", {
  d = rbind (c (0, 1, 2, 3), c (3, 2, 1, 0))
  expect_equal (fdm2id:::knn.indices (d, 2, drop.first = FALSE), matrix (c (1L, 4L, 2L, 3L),
                                                                         nrow = 2))
  expect_equal (fdm2id:::knn.indices (d, 2, drop.first = TRUE), matrix (c (2L, 3L, 3L, 2L),
                                                                        nrow = 2))
  # Per-row: the first row keeps its nearest, the second drops it.
  res = fdm2id:::knn.indices (d, 1, drop.first = c (FALSE, TRUE))
  expect_equal (as.vector (res), c (1L, 3L))
})

test_that ("selectfeatures() driven by Relief stays fast", {
  set.seed (0)
  n = 900
  x = do.call (rbind, lapply (1:3, function (i) matrix (stats::rnorm (n / 3 * 6, mean = i),
                                                        ncol = 6)))
  d = as.data.frame (x)
  labels = factor (rep (letters [1:3], each = n / 3))
  # Several seconds before the rewrite, well under one after it. Generous, so that a slow
  # machine does not fail the suite -- what is guarded against is a return to O(n^2 log n).
  expect_lt (system.time (suppressWarnings (
    selectfeatures (d, labels, algorithm = "ranking", unieval = "relief", uninb = 3)
  )) [["elapsed"]], 3)
})

# =========================================================================================
# Fourth audit, B8
# =========================================================================================

test_that ("CFS is built on a correlation, and is the same on every call", {
  data (iris)
  # r_cf used to be the Relief weight: not bounded the way Hall's formula assumes, possibly
  # negative, and estimated from a random sample -- so CFS answered differently every time.
  a = fdm2id:::fseval.cfs (iris [, -5], iris [, 5])
  b = fdm2id:::fseval.cfs (iris [, -5], iris [, 5])
  expect_equal (a, b)
  expect_true (is.finite (a))
  expect_gt (a, 0)
  # The merit of a single feature is its own correlation with the class, in [0, 1].
  one = fdm2id:::fseval.cfs (iris [, 3, drop = FALSE], iris [, 5])
  expect_gte (one, 0)
  expect_lte (one, 1)
  expect_equal (one, sqrt (unname (fdm2id:::fseval.inertiaratio (iris [, 3, drop = FALSE],
                                                                 iris [, 5]))))
  # Discriminant features must score above noise.
  set.seed (0)
  noise = data.frame (a = rnorm (150), b = rnorm (150))
  expect_gt (fdm2id:::fseval.cfs (iris [, 3:4], iris [, 5]),
             fdm2id:::fseval.cfs (noise, iris [, 5]))
})

test_that ("a constant variable separates nothing rather than dividing by zero", {
  data (iris)
  d = iris [, -5]
  d$flat = 1
  ratios = fdm2id:::fseval.inertiaratio (d, iris [, 5])
  expect_false (anyNA (ratios))
  expect_equal (unname (ratios ["flat"]), 0)
  expect_false (is.na (fdm2id:::fseval.cfs (d, iris [, 5])))
})

# =========================================================================================
# Fifth audit, D7 / D8 / D9
# =========================================================================================

test_that ("predict.selection picks its columns by name", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  model = suppressWarnings (FEATURESELECTION (d$train.x, d$train.y, uninb = 2,
                                              mainmethod = LDA))
  reference = predict (model, d$test.x)
  # Indexing by position silently used the wrong variables when the test set did not carry
  # its columns in the training order.
  expect_equal (as.character (predict (model, d$test.x [, 4:1])), as.character (reference))
  expect_equal (as.character (predict (model, d$test.x [, c (2, 4, 1, 3)])),
                as.character (reference))
  expect_error (predict (model, d$test.x [, 1:2]), "does not have the selected variable")
})

test_that ("the search algorithms do not leak their arguments into the learning method", {
  data (iris)
  # 'unieval', 'uninb' and 'unithreshold' used to fall into '...' and travel through the
  # evaluation criterion into the model being fitted.
  for (name in c ("fs.forward", "fs.backward", "fs.exhaustive", "fs.ranking"))
  {
    args = names (formals (get (name, envir = asNamespace ("fdm2id"))))
    for (a in c ("wrapmethod", "unieval", "uninb", "unithreshold"))
      expect_true (a %in% args, info = paste (name, a))
  }
  expect_error (suppressWarnings (
    selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "wrapper",
                    wrapmethod = LDA, nruns = 3)), NA)
})

test_that ("the wrapper's nruns reaches performance()", {
  # It was declared as 100 and never passed on, so performance() ran its own default of 10.
  expect_equal (formals (fdm2id:::fseval.wrapper)$nruns, 10)
  data (iris)
  fast = system.time (suppressWarnings (
    fdm2id:::fseval.wrapper (iris [, 3:4], iris [, 5], wrapmethod = LDA, nruns = 2)))
  slow = system.time (suppressWarnings (
    fdm2id:::fseval.wrapper (iris [, 3:4], iris [, 5], wrapmethod = LDA, nruns = 20)))
  expect_gt (slow [["elapsed"]], fast [["elapsed"]])
})

test_that ("a subset search does not return the top of the ranking", {
  data (iris)
  # A subset search drops Sepal.Length, which a ranking on Fisher's index places third, and
  # keeps Sepal.Width, which it places last: Sepal.Length is strongly correlated with both
  # petal variables and largely repeats them, where Sepal.Width is not and brings its own.
  selected = selectfeatures (iris [, -5], iris [, 5], algorithm = "forward",
                             multieval = "cfs")$features
  expect_setequal (selected, c ("Sepal.Width", "Petal.Length", "Petal.Width"))
})

# =========================================================================================
# Tenth audit: random forest importance, and drawing a selection
# =========================================================================================

test_that ("a random forest can rank the variables like any other univariate criterion", {
  skip_if_not_installed ("randomForest")
  data (iris)
  expect_true ("randomforest" %in% fdm2id:::fseval.univariate ())
  selection = selectfeatures (iris [, -5], iris [, 5], unieval = "randomforest", uninb = 2,
                              seed = 0)
  expect_s3_class (selection, "selection")
  expect_setequal (selection$features, c ("Petal.Length", "Petal.Width"))
  expect_named (selection$unieval, colnames (iris) [-5])
  # It goes through FEATURESELECTION like the others.
  expect_s3_class (FEATURESELECTION (iris [, -5], iris [, 5], unieval = "randomforest",
                                     uninb = 2, mainmethod = LDA, seed = 0), "selection")
  # And it works on a numeric target, where the importance is the decrease in MSE.
  data (trees)
  expect_length (selectfeatures (trees [, -3], trees [, 3], unieval = "randomforest",
                                 uninb = 1, seed = 0)$features, 1)
  # It scores variables one by one, so it refuses to score a subset.
  expect_message (out <- fdm2id:::fseval.randomforest (iris [, -5], iris [, 5],
                                                       vtype = "multivariate"), "univariate")
  expect_null (out)
})

test_that ("a selection made by ranking can be drawn", {
  skip_if_not_installed ("randomForest")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  selection = selectfeatures (iris [, -5], iris [, 5], unieval = "randomforest", uninb = 2,
                              seed = 0)
  before = graphics::par (c ("mar", "mfrow"))
  expect_error (plot (selection), NA)
  expect_equal (graphics::par (c ("mar", "mfrow")), before)
  expect_error (plot (selection, horiz = FALSE), NA)
  expect_error (plot (selection, legendpos = "topright"), NA)
  # A subset search scores subsets, not variables, so there is nothing to draw.
  subset = selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "fstat")
  expect_error (plot (subset), "no score per variable")
})

test_that ("a selection reports only the criteria it actually used", {
  data (iris)
  # Ranking with 'uninb' never calls the multivariate criterion, which was announced anyway.
  ranked = selectfeatures (iris [, -5], iris [, 5], unieval = "fisher", uninb = 2)
  expect_null (ranked$multivariate)
  expect_equal (ranked$univariate, "fisher")
  # Ranking that lets a multivariate criterion choose how many to keep does use it.
  sized = selectfeatures (iris [, -5], iris [, 5], unieval = "fisher", multieval = "fstat")
  expect_equal (sized$multivariate, "fstat")
})
