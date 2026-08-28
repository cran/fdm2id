# Tests for association.R
# arules/arulesViz are in Depends (always installed alongside fdm2id), so no
# skip_if_not_installed() guard is needed for them here.

test_that ("APRIORI smoke test: fit, predict, print, summary", {
  data (iris)
  d = arules::discretizeDF (iris,
      default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
  model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
  expect_s3_class (model, "apriori")
  pred = predict (model, d [, -5])
  expect_equal (length (pred), nrow (d))
  expect_output (print (model))
  expect_output (print (summary (model)))
})

test_that ("filter.rules() keeps only rules matching the requested side", {
  data ("Adult", package = "arules")
  r = arules::apriori (Adult, parameter = list (supp = .3, conf = .9), control = list (verbose = FALSE))
  right = filter.rules (r, right = "sex=")
  expect_true (length (right) <= length (r))
  expect_true (all (arules::`%pin%` (arules::rhs (right), "sex=")))
})

test_that ("general.rules() never returns more rules than it started with", {
  data ("Adult", package = "arules")
  r = arules::apriori (Adult, parameter = list (supp = .3, conf = .9), control = list (verbose = FALSE))
  reduced = general.rules (r)
  expect_true (length (reduced) <= length (r))
})

# --- Regression: predict.apriori()'s guard against an 'unmatched' label collision --------

test_that ("predict.apriori() rejects an 'unmatched' label that collides with an existing class", {
  data (iris)
  d = arules::discretizeDF (iris,
      default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
  model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
  expect_error (predict (model, d [, -5], unmatched = "setosa"))
})

# --- Regression: predict.apriori() shifted its labels ------------------------------------
# 'unmatched' was appended to the label vector only when length (unique (pred)) == n, so as
# soon as one consequent class was never predicted while some observations stayed unmatched,
# the vector was one entry short and those observations came out as NA.

test_that ("predict.apriori() never produces an NA label", {
  skip_if_not_installed ("arules")
  data (iris)
  d = arules::discretizeDF (iris, default = list (method = "interval", breaks = 3,
                                                  labels = c ("small", "medium", "large")))
  model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
  # Keep only the setosa and versicolor rules, then predict on setosa rows plus a handful of
  # rows no remaining rule covers: predictions are then {setosa, unmatched} out of three
  # possible codes -- exactly the case that used to mislabel.
  keep = arules::`%pin%` (arules::rhs (model$rules), "Class=setosa") |
         arules::`%pin%` (arules::rhs (model$rules), "Class=versicolor")
  model$rules = model$rules [keep]
  test = d [1:50, -5]
  weird = test [1:5, ]
  for (j in 1:4)
    weird [, j] = factor ("large", levels = levels (test [, j]))
  test = rbind (test, weird)
  pred = predict (model, test)
  expect_false (any (is.na (levels (pred))))
  expect_false (any (is.na (pred)))
  expect_true ("Unknown" %in% levels (pred))
})

test_that ("predict.apriori() reports a colliding 'unmatched' label in English", {
  skip_if_not_installed ("arules")
  data (iris)
  d = arules::discretizeDF (iris, default = list (method = "interval", breaks = 3,
                                                  labels = c ("small", "medium", "large")))
  model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
  expect_error (predict (model, d [, -5], unmatched = "setosa"), "collides")
})

# --- Fourth audit: predict.apriori() is reproducible ---------------------------------------

test_that ("predict.apriori() breaks ties without drawing lots", {
  require ("datasets")
  data (iris)
  d = discretizeDF (iris, default = list (method = "interval", breaks = 3,
                                          labels = c ("small", "medium", "large")))
  # A low confidence threshold leaves plenty of equally confident rules to choose between.
  model = APRIORI (d [, -5], d [, 5], supp = .02, conf = .5)
  set.seed (1)
  first = predict (model, d [, -5])
  set.seed (12345)
  second = predict (model, d [, -5])
  expect_equal (as.character (first), as.character (second))
  # And the choice does not depend on the state of the generator at all.
  expect_equal (as.character (predict (model, d [, -5])), as.character (first))
})

# --- Sixth audit: general.rules() keeps the subset matrix sparse ---------------------------

test_that ("general.rules() returns the same rules as the dense computation", {
  require ("arules")
  data (Adult)
  r = arules::apriori (Adult, parameter = list (supp = .3, conf = .8),
                       control = list (verbose = FALSE))
  dense = as.matrix (arules::is.subset (r@lhs, r@lhs) & arules::is.subset (r@rhs, r@rhs))
  diag (dense) = FALSE
  expected = r [!(colSums (dense, na.rm = TRUE) >= 1)]
  expect_equal (length (general.rules (r)), length (expected))
  expect_equal (arules::labels (general.rules (r)), arules::labels (expected))
})

# --- Tenth audit: arules' own parameters reach it ------------------------------------------

test_that ("APRIORI passes its extra arguments to arules", {
  require ("datasets")
  data (iris)
  d = discretizeDF (iris, default = list (method = "interval", breaks = 3,
                                          labels = c ("small", "medium", "large")))
  # '...' was declared and went nowhere, so 'maxlen' -- whose default of 10 truncates the
  # search and warns about it -- could not be reached.
  wide = APRIORI (d [, -5], d [, 5], supp = .05, conf = .7)
  short = APRIORI (d [, -5], d [, 5], supp = .05, conf = .7, maxlen = 3)
  expect_lt (length (short$rules), length (wide$rules))
  expect_length (predict (short, d [, -5]), nrow (d))
})
