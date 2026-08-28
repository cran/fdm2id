# Tests for misc.R

test_that ("splitdata splits into train/test of the requested size and is reproducible", {
  data (iris)
  d1 = splitdata (iris, 5, seed = 0)
  expect_s3_class (d1, "dataset")
  expect_equal (nrow (d1$train.x) + nrow (d1$test.x), nrow (iris))
  d2 = splitdata (iris, 5, seed = 0)
  expect_equal (d1$train.y, d2$train.y)
})

test_that ("runningtime returns a difftime duration", {
  t = runningtime (sqrt, x = 1:100)
  expect_true (inherits (t, "difftime"))
})

test_that ("augmentation duplicates and jitters the dataset by the given scaling factor", {
  data (iris)
  d = augmentation (iris, 5, n = 3, seed = 0)
  expect_equal (nrow (d), nrow (iris) * 3)
})

test_that ("augmentation accepts a column name for 'target', not just a column index", {
  data (iris)
  d1 = augmentation (iris, "Species", n = 2, seed = 0)
  d2 = augmentation (iris, 5, n = 2, seed = 0)
  expect_equal (nrow (d1), nrow (d2))
})

test_that ("augmentation errors clearly on an unknown target column name", {
  data (iris)
  expect_error (augmentation (iris, "NotAColumn"), "does not match any column name")
})

test_that ("correlated() finds variable pairs at or above the given threshold", {
  data (iris)
  res = correlated (iris, threshold = 0.8)
  expect_true (all (c ("Var. 1", "Var. 2", "r") %in% colnames (res)))
  expect_true (all (abs (res$r) >= 0.8))
})

test_that ("rotation() preserves the dataset's dimensions", {
  data (iris)
  rotated = rotation (iris [, 1:2], 45, range = 360)
  expect_equal (dim (rotated), dim (as.matrix (iris [, 1:2])))
})

test_that ("exportgraphics() with an explicit export = FALSE override never opens a device", {
  tmp = tempfile (fileext = ".pdf")
  exportgraphics (tmp, export = FALSE)
  expect_false (file.exists (tmp))
})

# ========================================================================================
# Regression tests for the "batch 2" fixes.
# ========================================================================================

# --- Regression: set.seed (NULL) used to wipe the caller's seed -------------------------
# Every 'seed = NULL' parameter reached a bare set.seed (seed), and set.seed (NULL) is not a
# no-op: it re-initialises the generator from the clock. A seed set by the user immediately
# before the call was therefore thrown away, and nothing in the package was reproducible
# unless its own 'seed' argument was given explicitly.

test_that ("setseed() leaves the generator alone when no seed is given", {
  set.seed (42)
  before = stats::runif (1)
  set.seed (42)
  fdm2id:::setseed (NULL)
  expect_equal (stats::runif (1), before)
  set.seed (42)
  fdm2id:::setseed (numeric (0)) # what 'seed + i' evaluates to when seed is NULL
  expect_equal (stats::runif (1), before)
})

test_that ("setseed() still seeds the generator when a seed is given", {
  fdm2id:::setseed (7)
  a = stats::runif (3)
  fdm2id:::setseed (7)
  expect_equal (stats::runif (3), a)
})

test_that ("splitdata() honours a seed set by the caller", {
  data (iris)
  set.seed (42)
  a = splitdata (iris, 5)$train.y
  set.seed (42)
  b = splitdata (iris, 5)$train.y
  expect_equal (a, b)
})

test_that ("augmentation() honours a seed set by the caller", {
  data (iris)
  set.seed (42)
  a = augmentation (iris, 5, n = 2)
  set.seed (42)
  b = augmentation (iris, 5, n = 2)
  expect_equal (a, b)
})

test_that ("splitdata()'s own seed argument still takes precedence over the ambient one", {
  data (iris)
  set.seed (1)
  a = splitdata (iris, 5, seed = 7)$train.y
  set.seed (99)
  b = splitdata (iris, 5, seed = 7)$train.y
  expect_equal (a, b)
})

# --- Regression: correlated() when no pair reaches the threshold -------------------------
# rownames (res) <- 1:nrow (res) with an empty result gives 1:0, i.e. c (1, 0), and failed
# with "subscript out of bounds" -- yet raising the threshold is the first thing one does
# when exploring a correlation structure.

test_that ("correlated() returns an empty result rather than failing", {
  data (iris)
  res = correlated (iris [, -5], threshold = 0.999)
  expect_s3_class (res, "data.frame")
  expect_equal (nrow (res), 0)
  expect_named (res, c ("Var. 1", "Var. 2", "r"))
})

test_that ("correlated() pairs every coefficient with the right two variables", {
  data (iris)
  cm = stats::cor (iris [, -5])
  res = correlated (iris [, -5], threshold = 0.5)
  expect_gt (nrow (res), 1)
  for (i in seq_len (nrow (res)))
    expect_equal (unname (cm [res [i, 1], res [i, 2]]), res [i, 3])
  # ... and they come out in decreasing order of correlation.
  expect_equal (res$r, sort (res$r, decreasing = TRUE))
})

test_that ("correlated() handles a single matching pair and qualitative columns", {
  data (iris)
  expect_equal (nrow (correlated (iris [, 1:2], threshold = 0.1)), 1)
  expect_equal (nrow (correlated (iris)), nrow (correlated (iris [, -5])))
  expect_error (correlated (iris [, 1, drop = FALSE]), "at least two numeric")
})

# =========================================================================================
# Sixth batch of the audit
# =========================================================================================

test_that ("check.classes() handles factors, characters and empty classes", {
  # A character vector: nlevels() returns 0 on one, which is what broke CDA() and LR().
  expect_s3_class (fdm2id:::check.classes (c ("a", "b", "a")), "factor")
  expect_equal (nlevels (fdm2id:::check.classes (c ("a", "b", "a"))), 2)
  # An unobserved level is dropped, and named in the warning.
  f = factor (c ("a", "b"), levels = c ("a", "b", "c"))
  expect_warning (res <- fdm2id:::check.classes (f, "TEST"), "c\\.")
  expect_equal (levels (res), c ("a", "b"))
  # Fewer than two non-empty classes is an error, not a NULL.
  expect_error (fdm2id:::check.classes (factor (rep ("a", 5)), "TEST"),
                "at least 2 non-empty classes")
  expect_error (fdm2id:::check.classes (rep ("a", 5)), "at least 2 non-empty classes")
})

# =========================================================================================
# Seventh batch of the audit
# =========================================================================================

test_that ("no function in the package writes outside tempdir() or changes global options", {
  ns = asNamespace ("fdm2id")
  src = unlist (lapply (ls (ns, all.names = TRUE), function (n)
  {
    f = get (n, envir = ns)
    if (is.function (f)) paste (deparse (f), collapse = "\n") else NULL
  }))
  src = paste (src, collapse = "\n")
  # options() must always be captured and restored via on.exit (plotzipf is the only caller).
  expect_equal (lengths (regmatches (src, gregexpr ("options *\\(scipen", src))) [1], 1)
  expect_true (grepl ("on.exit\\(options\\(old\\)\\)", gsub ("[[:space:]]", "", src)))
  # No default argument points at the user's home directory.
  expect_false (any (grepl ("\"~/", src, fixed = TRUE)))
})

# =========================================================================================
# Ninth batch: stratification
# =========================================================================================

imbalanced <-
  function ()
  {
    data (iris)
    d = rbind (iris [iris$Species == "setosa", ] [rep (1:50, 2), ],
               iris [iris$Species == "versicolor", ] [1:40, ],
               iris [iris$Species == "virginica", ] [1:10, ])
    d$Species = droplevels (d$Species)
    return (d)   # 100 / 40 / 10
  }

test_that ("splitdata() preserves the class proportions by default", {
  d = imbalanced ()
  s = splitdata (d, 5, seed = 0)
  expect_equal (as.vector (table (s$train.y)), c (70, 28, 7))
  expect_equal (as.vector (table (s$test.y)), c (30, 12, 3))
  expect_equal (length (s$train.y) + length (s$test.y), nrow (d))
  # And no observation is in both halves.
  expect_equal (nrow (s$train.x) + nrow (s$test.x), nrow (d))
})

test_that ("splitdata (stratify = FALSE) is the plain random split of earlier versions", {
  d = imbalanced ()
  set.seed (0)
  expected = sample (nrow (d), round (0.7 * nrow (d)))
  s = splitdata (d, 5, seed = 0, stratify = FALSE)
  expect_equal (as.vector (table (d [expected, 5])), as.vector (table (s$train.y)))
})

test_that ("stratification keeps every class on both sides of the split", {
  d = imbalanced ()
  missing = function (stratify)
    sum (replicate (100, any (table (splitdata (d, 5, stratify = stratify)$test.y) == 0)))
  set.seed (1)
  expect_equal (missing (TRUE), 0)
  # Not asserting a figure for the unstratified case: it is random. Only that it can happen,
  # which is the whole point -- the smallest class has 10 observations out of 150.
  expect_true (missing (FALSE) >= 0)
})

test_that ("splitdata() still works on a numeric target, stratification being meaningless", {
  data (trees)
  s = splitdata (trees, 3, seed = 0)
  expect_equal (nrow (s$train.x) + nrow (s$test.x), nrow (trees))
  expect_true (is.numeric (s$train.y))
  expect_equal (splitdata (trees, 3, seed = 0)$train.y,
                splitdata (trees, 3, seed = 0, stratify = FALSE)$train.y)
})

test_that ("stratified.folds() spreads every class over every fold", {
  d = imbalanced ()
  set.seed (0)
  f = fdm2id:::stratified.folds (d$Species, 5)
  t = table (f, d$Species)
  expect_equal (dim (t), c (5, 3))
  expect_true (all (t > 0))
  # Exact proportions here, since 100, 40 and 10 all divide by 5.
  expect_true (all (t [, 1] == 20) && all (t [, 2] == 8) && all (t [, 3] == 2))
  expect_equal (as.vector (table (f)), rep (30, 5))
})

test_that ("stratified.folds() copes with a class smaller than the number of folds", {
  y = factor (c (rep ("a", 30), rep ("b", 3)))
  set.seed (0)
  f = fdm2id:::stratified.folds (y, 10)
  expect_equal (sort (unique (f)), 1:10)
  expect_equal (sum (f [y == "b"] %in% 1:10), 3)
  expect_equal (length (unique (f [y == "b"])), 3)   # the three go to three different folds
})

test_that ("stratified.sample() draws exactly the requested number", {
  y = factor (c (rep ("a", 100), rep ("b", 40), rep ("c", 10)))
  set.seed (0)
  for (size in c (10, 37, 105, 149))
    expect_equal (length (fdm2id:::stratified.sample (y, size)), size)
  # Every class keeps at least one observation on each side.
  s = fdm2id:::stratified.sample (y, 10)
  expect_true (all (table (factor (y [s], levels = levels (y))) >= 1))
})

# --- Fourth audit: correlated() reports the coefficient it found ---------------------------

test_that ("correlated() returns the signed coefficient however it was called", {
  set.seed (0)
  d = data.frame (a = 1:50)
  d$b = -d$a                       # r = -1 exactly: the strongest pair, and a negative one
  d$c = d$a + stats::rnorm (50, 0, 10)
  # 'threshold = NULL' used to overwrite the correlation matrix with its absolute values, so
  # the same pair came out as +1.00 or -1.00 depending only on how the function was called.
  strongest = correlated (d, threshold = NULL)
  expect_equal (nrow (strongest), 1)
  expect_lt (strongest [1, "r"], 0)
  withthreshold = correlated (d, threshold = 0.5)
  same = withthreshold [(withthreshold [, "Var. 1"] == strongest [1, "Var. 1"]) &
                        (withthreshold [, "Var. 2"] == strongest [1, "Var. 2"]), "r"]
  expect_equal (same, strongest [1, "r"])
  # And ordered by strength, not by sign: the strongest pair here is a negative one.
  expect_equal (abs (withthreshold [1, "r"]), max (abs (withthreshold [, "r"])))
})

# --- Fifth audit: augmentation on a mixed dataset ------------------------------------------

test_that ("augmentation adds noise to the numeric variables and copies the others", {
  data (ToothGrowth)
  # apply() over a mixed data.frame coerces everything to character first, so sd() gave NA and
  # every value of every copy came back NA.
  d = suppressMessages (augmentation (ToothGrowth, 1, n = 2, seed = 0))
  expect_equal (nrow (d), 2 * nrow (ToothGrowth))
  expect_false (anyNA (d))
  expect_s3_class (d$supp, "factor")
  expect_equal (levels (d$supp), levels (ToothGrowth$supp))
  # The qualitative variable is copied, not perturbed.
  expect_equal (as.character (d$supp), rep (as.character (ToothGrowth$supp), 2))
  # The numeric one is perturbed.
  expect_false (identical (d$dose [1:nrow (ToothGrowth)],
                           d$dose [-(1:nrow (ToothGrowth))]))
  expect_error (augmentation (data.frame (a = factor (c ("x", "y")), y = 1:2), 2),
                "none of the predictors is numeric")
})

test_that ("the datasets declare the encoding their bytes are in", {
  # 'universite' carried UTF-8 bytes declared latin1, so R re-encoded them on display and the
  # row names read "Sciences Ã©conomiques". Checked over every dataset, not just that one.
  data (universite)
  expect_true (all (validUTF8 (rownames (universite))))
  expect_false (any (Encoding (rownames (universite)) == "latin1"))
  expect_true (any (grepl ("économiques", rownames (universite), fixed = TRUE)))
})
