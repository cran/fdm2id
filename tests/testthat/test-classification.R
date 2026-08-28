# Tests for classification.R
#
# These are deliberately modest smoke tests (does it run, does it return the right shape),
# plus a handful of regression tests for bugs found and fixed during the Phase 3 review
# sessions -- they are the highest-value tests to have, since they would have caught these
# bugs in seconds instead of a full R CMD check.

test_that ("NB / LDA smoke tests: fit + predict return the expected shape", {
  skip_if_not_installed ("e1071")
  skip_if_not_installed ("MASS")
  data (iris)
  for (ctor in list (NB, LDA))
  {
    model = ctor (iris [, -5], iris [, 5])
    expect_s3_class (model, "model")
    pred = predict (model, iris [, -5])
    expect_equal (length (pred), nrow (iris))
    expect_true (is.factor (pred))
  }
})

# --- Regression: xgboost >= 2.1 interface (x/y/learning_rate/verbosity) ------------------
# GRADIENTBOOSTING used to recode labels to 0-based numeric and pass objective =
# "multi:softprob" by hand; recent xgboost releases reject that and require a factor 'y'
# with the objective inferred automatically. See NEWS.md for details.

test_that ("GRADIENTBOOSTING works with the xgboost >= 2.1 interface (multiclass)", {
  skip_if_not_installed ("xgboost")
  data (iris)
  model = GRADIENTBOOSTING (iris [, -5], iris [, 5], ntree = 10)
  pred = predict (model, iris [, -5])
  expect_true (is.factor (pred))
  expect_equal (levels (pred), levels (iris [, 5]))
  fuzzy = predict (model, iris [, -5], fuzzy = TRUE)
  expect_true (is.matrix (fuzzy))
  expect_equal (ncol (fuzzy), nlevels (iris [, 5]))
})

test_that ("GRADIENTBOOSTING works with the xgboost >= 2.1 interface (binary)", {
  skip_if_not_installed ("xgboost")
  d = iris
  levels (d [, 5]) = c ("+", "+", "-") # two-class problem
  model = GRADIENTBOOSTING (d [, -5], d [, 5], ntree = 10)
  pred = predict (model, d [, -5])
  expect_true (is.factor (pred))
  fuzzy = predict (model, d [, -5], fuzzy = TRUE)
  expect_true (is.matrix (fuzzy))
  expect_equal (ncol (fuzzy), 2)
  expect_equal (colnames (fuzzy), levels (d [, 5]))
})

# --- Regression: RANDOMFOREST seed + methodparameters/graph absorption -------------------

test_that ("RANDOMFOREST seed gives reproducible predictions", {
  skip_if_not_installed ("randomForest")
  data (iris)
  m1 = RANDOMFOREST (iris [, -5], iris [, 5], seed = 0)
  m2 = RANDOMFOREST (iris [, -5], iris [, 5], seed = 0)
  expect_equal (predict (m1, iris [, -5]), predict (m2, iris [, -5]))
})

test_that ("RANDOMFOREST absorbs graph/methodparameters instead of forwarding them to randomForest::randomForest()", {
  skip_if_not_installed ("randomForest")
  data (iris)
  # This is exactly the call shape performance()/protocol.fitpredict() uses internally.
  expect_error (
    RANDOMFOREST (iris [, -5], iris [, 5], graph = FALSE, methodparameters = NULL, seed = 0),
    NA
  )
})

# --- Regression: protocol.holdout() single-predictor is.vector() copy-paste bug ----------
# The second is.vector() guard used to re-test train.x instead of test.x, so a single-column
# test.x collapsed to a bare vector by [-s, ] indexing was never re-wrapped into a data.frame.

test_that ("performance() with protocol = 'holdout' works on a single-predictor dataset", {
  skip_if_not_installed ("e1071")
  data (iris)
  expect_error (
    performance (NB, iris [, 1], iris [, 5], protocol = "holdout", seed = 0),
    NA
  )
})

test_that ("performance() with protocol = 'holdout' still works on a multi-predictor dataset", {
  skip_if_not_installed ("e1071")
  data (iris)
  expect_error (
    performance (NB, iris [, -5], iris [, 5], protocol = "holdout", seed = 0),
    NA
  )
})

# --- Regression: evaluation() dispatch table (get(paste(...)) -> named list) -------------

test_that ("evaluation() dispatch table still resolves eval.* functions correctly", {
  skip_if_not_installed ("e1071")
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  model = NB (d$train.x, d$train.y)
  pred = predict (model, d$test.x)
  res = evaluation (pred, d$test.y, eval = c ("accuracy", "kappa"))
  expect_named (res, c ("accuracy", "kappa"))
  expect_true (all (res >= -1 & res <= 1))
})

test_that ("evaluation() also dispatches to regression criteria defined in regression.R", {
  data (trees)
  d = splitdata (trees, 3, seed = 0)
  model = LINREG (d$train.x, d$train.y)
  pred = predict (model, d$test.x)
  res = evaluation (pred, d$test.y, eval = c ("r2", "msep"))
  expect_named (res, c ("r2", "msep"))
})

test_that ("evaluation() gives a clear error for an unknown criterion", {
  data (iris)
  expect_error (evaluation (iris [, 5], iris [, 5], eval = "not_a_real_criterion"),
               "unknown evaluation criterion")
})

# --- Regression: predict.model() dispatch table (if/else chain -> named list) ------------

test_that ("predict.model() dispatch table covers CART/LDA/QDA/NB alike", {
  skip_if_not_installed ("e1071")
  skip_if_not_installed ("MASS")
  skip_if_not_installed ("rpart")
  data (iris)
  for (ctor in list (NB, LDA, QDA, CART))
  {
    model = ctor (iris [, -5], iris [, 5])
    pred = predict (model, iris [, -5])
    expect_equal (length (pred), nrow (iris))
  }
})

test_that ("predict.model() falls back to generic stats::predict() for unrecognised methods", {
  data (trees)
  # Build a "model"-class object whose $method is NOT in predictmodel.functions()'s table,
  # to directly exercise the fallback branch (rather than one of the ~13 dispatched methods).
  fake = list (model = stats::lm (Volume ~ Girth, trees), method = "SOME_UNKNOWN_METHOD")
  class (fake) = "model"
  expect_error (predict (fake, trees), NA)
})

# ========================================================================================
# Regression tests for the "batch 1" fixes (evaluation measures returning wrong values).
#
# The three bugs below all returned an object of the *right shape* with the *wrong value*,
# which is exactly what a smoke test cannot catch. Each test therefore compares against a
# figure computed by hand, or against an invariant that must hold whatever the method.
# ========================================================================================

# --- Regression: evaluation.precision() / evaluation.recall() argument roles -------------
# Both functions used to divide by the wrong margin of table (gt, predictions): precision
# returned the recall and vice versa. Reading the code alone could not settle it, because
# calling them ground-truth-first (the order the confusion() example used to show) gave the
# right numbers -- but evaluation() and performance() bind these arguments by NAME, so
# through the high-level API the inversion was unavoidable.

test_that ("evaluation.precision() and evaluation.recall() match hand-computed values", {
  # 4 actual positives, 3 predicted positives, 2 true positives.
  gt   = factor (c ("+", "+", "+", "+", "-", "-", "-", "-", "-", "-"))
  pred = factor (c ("+", "+", "-", "-", "+", "-", "-", "-", "-", "-"))
  # levels(gt) is c("+", "-"), so the default positive class is "+".
  expect_equal (evaluation.precision (pred, gt), 2 / 3) # TP / predicted positives
  expect_equal (evaluation.recall    (pred, gt), 2 / 4) # TP / actual positives
  # And with the other class as the positive one: 5 TN, 6 actual "-", 7 predicted "-".
  expect_equal (evaluation.precision (pred, gt, positive = "-"), 5 / 7)
  expect_equal (evaluation.recall    (pred, gt, positive = "-"), 5 / 6)
})

test_that ("precision and recall keep their roles through evaluation() and performance()", {
  gt   = factor (c ("+", "+", "+", "+", "-", "-", "-", "-", "-", "-"))
  pred = factor (c ("+", "+", "-", "-", "+", "-", "-", "-", "-", "-"))
  expect_equal (unname (evaluation (pred, gt, eval = "precision")), 2 / 3)
  expect_equal (unname (evaluation (pred, gt, eval = "recall")),    2 / 4)
  # The F-measure is symmetric in precision and recall, which is precisely why the
  # inversion stayed invisible for so long -- check it explicitly all the same.
  expect_equal (unname (evaluation (pred, gt, eval = "fmeasure")),
                2 * (2 / 3) * (2 / 4) / ((2 / 3) + (2 / 4)))
})

test_that ("evaluation.recall() names itself in its error message", {
  data (iris)
  # Three classes is no longer an error (ninth batch): it is macro-averaged. Asking for the
  # two-class reading explicitly still is, and the message still names the function.
  expect_error (evaluation.recall (iris [, 5], iris [, 5], average = "binary"),
                "evaluation.recall")
})

# --- Regression: accuracy on a non-square contingency table ------------------------------
# eval.accuracy() used sum (diag (table (p, l))): table() only keeps observed values, so as
# soon as predictions and ground truth had different supports the diagonal read unrelated
# cells and returned a plausible but wrong figure -- up to 1 for a completely wrong model.

test_that ("evaluation.accuracy() is correct when predictions and ground truth have different supports", {
  gt   = factor (c ("a", "a", "b", "b"), levels = c ("a", "b", "c"))
  pred = factor (c ("b", "b", "c", "c"), levels = c ("a", "b", "c"))
  expect_equal (evaluation.accuracy (pred, gt), 0) # not a single correct prediction
})

test_that ("evaluation.accuracy() is correct when a class is never predicted", {
  gt   = factor (c ("a", "a", "b", "b", "c", "c"))
  pred = factor (c ("a", "a", "b", "b", "b", "b"), levels = c ("a", "b", "c"))
  expect_equal (evaluation.accuracy (pred, gt), 4 / 6)
})

test_that ("evaluation.accuracy() still agrees with the naive count in the ordinary case", {
  data (iris)
  set.seed (0)
  pred = sample (iris [, 5])
  expect_equal (evaluation.accuracy (pred, iris [, 5]),
                sum (pred == iris [, 5]) / nrow (iris))
})

# --- Regression: fuzzy predictions of the two-class LR / MLP models ----------------------
# nnet::multinom (type = "probs") and nnet::nnet (type = "raw") both return the probability
# of the *second* level only for a two-class problem. The wrapper did cbind (res, 1 - res)
# and then named the columns after the levels, so the columns ended up swapped: argmax of
# the fuzzy prediction contradicted the hard prediction on *every* observation.

test_that ("argmax of the fuzzy prediction agrees with the hard prediction (two classes)", {
  skip_if_not_installed ("e1071")
  skip_if_not_installed ("MASS")
  skip_if_not_installed ("nnet")
  data (iris)
  d = iris
  levels (d [, 5]) = c ("A", "A", "B")
  models = list (LR  = LR  (d [, -5], d [, 5]),
                 MLP = MLP (d [, -5], d [, 5], hidden = 3, decay = .1),
                 NB  = NB  (d [, -5], d [, 5]),
                 LDA = LDA (d [, -5], d [, 5]))
  for (name in names (models))
  {
    fuzzy = predict (models [[name]], d [, -5], fuzzy = TRUE)
    hard  = as.character (predict (models [[name]], d [, -5]))
    expect_equal (colnames (fuzzy), levels (d [, 5]), info = name)
    expect_equal (colnames (fuzzy) [apply (fuzzy, 1, which.max)], hard, info = name)
  }
})

test_that ("argmax of the fuzzy prediction agrees with the hard prediction (three classes)", {
  skip_if_not_installed ("e1071")
  skip_if_not_installed ("nnet")
  data (iris)
  models = list (LR  = LR  (iris [, -5], iris [, 5]),
                 MLP = MLP (iris [, -5], iris [, 5], hidden = 4, decay = .1),
                 NB  = NB  (iris [, -5], iris [, 5]))
  for (name in names (models))
  {
    fuzzy = predict (models [[name]], iris [, -5], fuzzy = TRUE)
    hard  = as.character (predict (models [[name]], iris [, -5]))
    expect_equal (colnames (fuzzy) [apply (fuzzy, 1, which.max)], hard, info = name)
  }
})

test_that ("the two-class LR probabilities match nnet::multinom's own output", {
  skip_if_not_installed ("nnet")
  data (iris)
  d = iris
  levels (d [, 5]) = c ("A", "A", "B")
  reference = stats::predict (nnet::multinom (Class ~ .,
                                              cbind.data.frame (d [, -5], Class = d [, 5]),
                                              trace = FALSE),
                              d [, -5], type = "probs")
  fuzzy = predict (LR (d [, -5], d [, 5]), d [, -5], fuzzy = TRUE)
  # 'reference' is P(second level) = P("B").
  expect_equal (unname (fuzzy [, "B"]), unname (reference))
})

# --- Regression: confusion() orientation -------------------------------------------------
# The matrix is documented (and computed) with the true labels on the rows; the axes of the
# plot used to be labelled the other way round, and the example passed the ground truth
# first, which silently transposed the result.

test_that ("confusion() puts the true labels on the rows, as documented", {
  gt   = factor (c ("+", "+", "+", "+", "-", "-", "-", "-", "-", "-"))
  pred = factor (c ("+", "+", "-", "-", "+", "-", "-", "-", "-", "-"))
  conf = confusion (pred, gt, norm = FALSE, graph = FALSE)
  expect_equal (names (dimnames (conf)), c ("True labels", "Predicted labels"))
  expect_equal (conf ["+", "+"], 2) # true positives
  expect_equal (conf ["-", "+"], 1) # false positives: actually "-", predicted "+"
  expect_equal (conf ["+", "-"], 2) # false negatives
})

test_that ("confusion() stays square when a class is never predicted", {
  gt   = factor (c ("a", "a", "b", "b", "c", "c"))
  pred = factor (c ("a", "a", "b", "b", "b", "b"), levels = c ("a", "b", "c"))
  conf = confusion (pred, gt, norm = FALSE, graph = FALSE)
  expect_equal (dim (conf), c (3L, 3L))
  expect_equal (rownames (conf), colnames (conf))
})

test_that ("confusion() normalisation does not produce NaN for an unobserved class", {
  gt   = factor (c ("a", "a", "b", "b"), levels = c ("a", "b", "c"))
  pred = factor (c ("a", "a", "b", "c"), levels = c ("a", "b", "c"))
  conf = confusion (pred, gt, norm = TRUE, graph = FALSE)
  expect_false (any (is.nan (conf)))
  expect_equal (unname (conf ["a", "a"]), 1)
})

# --- Regression: MLP()'s default 'hidden' was collapsed by ifelse() ----------------------
# ifelse() is vectorised: it returns only the first element of the branch it selects, so the
# default hidden = ifelse (is.vector (train), 2:(1 + k), 2:(ncol + k)) always evaluated to
# the single value 2, and the cross-validated choice of the hidden layer size documented for
# this argument never took place.

test_that ("MLP()'s default hidden layer size is a grid, not a single value", {
  data (iris)
  grid = eval (formals (MLP)$hidden, list (train = iris [, -5], labels = iris [, 5]))
  expect_equal (grid, 2:(ncol (iris [, -5]) + nlevels (iris [, 5])))
  expect_gt (length (grid), 1)
  # ... and the single-predictor branch too.
  grid = eval (formals (MLP)$hidden, list (train = iris [, 1], labels = iris [, 5]))
  expect_equal (grid, 2:(1 + nlevels (iris [, 5])))
  expect_gt (length (grid), 1)
})

test_that ("MLPREG()'s default size is unchanged by the ifelse() -> if/else rewrite", {
  data (trees)
  expect_equal (eval (formals (MLPREG)$size, list (x = trees [, -3])), 2:ncol (trees [, -3]))
  expect_equal (eval (formals (MLPREG)$size, list (x = trees [, 1])), 2)
})

test_that ("performance() is reproducible under a caller-set seed", {
  skip_if_not_installed ("e1071")
  data (iris)
  set.seed (42)
  a = performance (NB, iris [, -5], iris [, 5], nruns = 3)
  set.seed (42)
  b = performance (NB, iris [, -5], iris [, 5], nruns = 3)
  expect_equal (a, b)
})

# ========================================================================================
# Regression tests for the "batch 3" fixes (things that crashed, or did nothing silently).
# ========================================================================================

# --- Regression: LDA/QDA with a single predictor -----------------------------------------
# The is.vector() guard was immediately undone by train [, variables], which drops the
# dimension when a single column survives the zero-variance filter. MASS then rejected the
# bare vector with "'x' is not a matrix", making univariate discriminant analysis impossible
# -- and breaking selectfeatures (..., multieval = "wrapper"), which necessarily evaluates
# one-variable subsets.

test_that ("LDA() and QDA() work with a single predictor", {
  skip_if_not_installed ("MASS")
  data (iris)
  for (ctor in list (LDA, QDA))
  {
    model = ctor (iris [, 1, drop = FALSE], iris [, 5])
    expect_s3_class (model, "model")
    expect_equal (length (predict (model, iris [, 1, drop = FALSE])), nrow (iris))
    # ... and given a bare vector, which the is.vector() guard is meant to handle.
    model = ctor (iris [, 1], iris [, 5])
    expect_equal (length (predict (model, iris [, 1])), nrow (iris))
  }
})

test_that ("LDA() predictions on a single predictor agree with MASS::lda() called directly", {
  skip_if_not_installed ("MASS")
  data (iris)
  reference = stats::predict (MASS::lda (x = iris [, 1, drop = FALSE], grouping = iris [, 5]),
                              iris [, 1, drop = FALSE])$class
  expect_equal (predict (LDA (iris [, 1, drop = FALSE], iris [, 5]), iris [, 1, drop = FALSE]),
                reference)
})

test_that ("LDA() drops constant predictors and still predicts", {
  skip_if_not_installed ("MASS")
  data (iris)
  d = iris [, 1:3]
  d [, 2] = 1 # zero variance
  model = LDA (d, iris [, 5])
  expect_equal (sum (model$variables), 2)
  expect_equal (length (predict (model, d)), nrow (iris))
})

test_that ("LDA() reports non-numeric predictors clearly", {
  skip_if_not_installed ("MASS")
  data (iris)
  d = iris [, -5]
  d$flag = factor (rep (c ("a", "b"), 75))
  expect_error (LDA (d, iris [, 5]), "must be numeric")
  expect_error (LDA (data.frame (a = rep (1, 10), b = rep (2, 10)),
                     factor (rep (c ("x", "y"), 5))), "constant")
})

# --- Regression: evaluation() used to answer NULL for precision/recall beyond two classes --
# It then raised an explicit error (third batch); since the ninth batch it computes the
# macro-averaged value, which is what one wanted in the first place.

test_that ("evaluation() computes precision and recall beyond two classes", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pred = predict (NB (d$train.x, d$train.y), d$test.x)
  res = evaluation (pred, d$test.y, eval = c ("precision", "recall"))
  expect_named (res, c ("precision", "recall"))
  expect_equal (unname (res),
                c (evaluation.precision (pred, d$test.y, average = "macro"),
                   evaluation.recall (pred, d$test.y, average = "macro")))
  # Perfect predictions score 1 whatever the number of classes.
  expect_equal (unname (evaluation (iris [, 5], iris [, 5], eval = "precision")), 1)
})

test_that ("evaluation() still returns precision and recall on two classes", {
  d = iris
  levels (d [, 5]) = c ("A", "A", "B")
  res = evaluation (d [, 5], d [, 5], eval = c ("precision", "recall"))
  expect_named (res, c ("precision", "recall"))
  expect_equal (unname (res), c (1, 1))
})

# --- Regression: cartdepth() on a tree reduced to its root -------------------------------
# ceiling (log2 (max node)) - 1 agreed with floor (log2 (max node)) everywhere except on a
# single-node tree, where it returned -1.

test_that ("cartdepth() returns 0 for a tree that could not be split", {
  skip_if_not_installed ("rpart")
  model = CART (data.frame (X = c (1, 1, 1, 1)), factor (c ("a", "b", "a", "b")))
  expect_equal (cartnodes (model), 1)
  expect_equal (cartdepth (model), 0)
})

test_that ("cartdepth() matches the depth read off rpart's node numbering", {
  skip_if_not_installed ("rpart")
  data (iris)
  model = CART (iris [, -5], iris [, 5])
  nodes = as.numeric (rownames (model$model$frame))
  expect_equal (cartdepth (model), floor (log2 (max (nodes))))
})

# --- Regression: CART() cross-validated cp with xval = nrow (d) ---------------------------
# One tree fitted per observation, purely to choose cp: quadratic, 9.8 s on 5000 rows.

test_that ("CART() exposes xval and defaults to 10 folds", {
  skip_if_not_installed ("rpart")
  data (iris)
  expect_equal (eval (formals (CART)$xval), 10)
  # The old behaviour is still reachable, and gives the same tree on iris.
  a = CART (iris [, -5], iris [, 5])
  b = CART (iris [, -5], iris [, 5], xval = nrow (iris))
  expect_equal (cartleafs (a), cartleafs (b))
  expect_equal (predict (a, iris [, -5]), predict (b, iris [, -5]))
})

test_that ("CART() with an explicit cp ignores xval entirely", {
  skip_if_not_installed ("rpart")
  data (iris)
  expect_error (CART (iris [, -5], iris [, 5], cp = 0.01, xval = 0), NA)
})

# =========================================================================================
# Sixth batch of the audit: API contract and design residues
# =========================================================================================

# --- Regression: ROC / cost curves drawn from hard labels ---------------------------------
# roc.curves() used to be handed the predicted *labels*, which collapse a ROC curve to a
# single operating point joined to (0,0) and (1,1). Both readings are now available, and
# performance() defaults to the one built on class-membership scores.

two.class.data <-
  function ()
  {
    data (iris)
    d = iris [iris [, 5] != "versicolor", ]
    d [, 5] = factor (as.character (d [, 5]), levels = c ("virginica", "setosa"))
    return (d)
  }

test_that ("curve.scores() reads hard labels and scores, and refuses the impossible one", {
  gt = factor (c ("a", "a", "b", "b"))
  # Hard labels -> 0/1 indicator of the positive class.
  expect_equal (fdm2id:::curve.scores (factor (c ("a", "b", "b", "b")), gt, "a", "hard") [[1]],
                c (1, 0, 0, 0))
  # A numeric vector is taken as scores.
  s = c (.9, .8, .2, .1)
  expect_equal (fdm2id:::curve.scores (s, gt, "a", "fuzzy") [[1]], s)
  # A matrix with one column per class: the positive column is picked by name, not position.
  m = cbind (a = s, b = 1 - s)
  expect_equal (fdm2id:::curve.scores (m, gt, "b", "fuzzy") [[1]], 1 - s)
  expect_equal (fdm2id:::curve.scores (m, gt, "a", "fuzzy") [[1]], s)
  # Asking for a continuous curve from labels cannot work, and says so.
  expect_error (fdm2id:::curve.scores (factor (c ("a", "b", "b", "b")), gt, "a", "fuzzy"),
                "hard class labels")
})

test_that ("curve.scores (type = 'auto') picks the right reading", {
  gt = factor (c ("a", "a", "b", "b"))
  expect_equal (fdm2id:::curve.scores (factor (c ("a", "b", "b", "b")), gt, "a") [[1]],
                c (1, 0, 0, 0))
  expect_equal (fdm2id:::curve.scores (c (.9, .8, .2, .1), gt, "a") [[1]], c (.9, .8, .2, .1))
})

test_that ("roc.curves() accepts both modes and validates its arguments", {
  skip_if_not_installed ("ROCR")
  d = two.class.data ()
  model = NB (d [, -5], d [, 5])
  hard = predict (model, d [, -5])
  fuzzy = predict (model, d [, -5], fuzzy = TRUE)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (roc.curves (hard, d [, 5], type = "hard"), NA)
  expect_error (roc.curves (fuzzy, d [, 5], type = "fuzzy"), NA)
  expect_error (roc.curves (fuzzy, d [, 5], positive = "setosa"), NA)
  # A three-class problem has no ROC curve.
  expect_error (roc.curves (iris [, 5], iris [, 5]), "two-class")
  # A positive class that does not exist is caught before ROCR sees it.
  expect_error (roc.curves (fuzzy, d [, 5], positive = "versicolor"), "not one of the class")
})

test_that ("cost.curves() accepts both modes", {
  skip_if_not_installed ("ROCR")
  d = two.class.data ()
  model = NB (d [, -5], d [, 5])
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (cost.curves (predict (model, d [, -5]), d [, 5], type = "hard"), NA)
  expect_error (cost.curves (predict (model, d [, -5], fuzzy = TRUE), d [, 5],
                             type = "fuzzy"), NA)
})

test_that ("the positive class drives the orientation of the curve", {
  skip_if_not_installed ("ROCR")
  gt = factor (c ("neg", "neg", "pos", "pos"), levels = c ("neg", "pos"))
  s = c (.1, .2, .8, .9)   # a perfect ranking for "pos"
  auc = function (positive, scores)
  {
    sc = fdm2id:::curve.scores (scores, gt, positive, "fuzzy") [[1]]
    p = ROCR::prediction (sc, as.numeric (gt == positive))
    ROCR::performance (p, "auc")@y.values [[1]]
  }
  expect_equal (auc ("pos", s), 1)
  expect_equal (auc ("neg", s), 0)   # the mirrored curve the old code could draw
})

# --- Regression: the five protocol.* functions had five different signatures ---------------

test_that ("every protocol.* function accepts the same argument list", {
  common = c ("methods", "train.x", "train.y", "test.x", "test.y", "train.size",
              "methodparameters", "nruns", "nfolds", "seed", "fuzzy")
  for (name in names (fdm2id:::protocol.functions ()))
  {
    f = fdm2id:::protocol.functions () [[name]]
    expect_true (all (common %in% names (formals (f))),
                 info = paste ("protocol", name, "is missing:",
                               paste (setdiff (common, names (formals (f))), collapse = ", ")))
  }
})

test_that ("performance() runs each protocol without an 'unused argument' error", {
  data (iris)
  d = iris [iris [, 5] != "versicolor", ]
  d [, 5] = factor (as.character (d [, 5]))
  for (p in c ("train", "holdout", "crossvalidation", "bootstrap"))
    expect_error (performance (NB, d [, -5], d [, 5], protocol = p, nruns = 2, nfolds = 3,
                               seed = 0), NA,
                  info = paste ("protocol =", p))
})

test_that ("performance() draws ROC curves from hard labels unless fuzzy is asked for", {
  skip_if_not_installed ("ROCR")
  data (iris)
  d = iris [iris [, 5] != "versicolor", ]
  d [, 5] = factor (as.character (d [, 5]))
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_equal (eval (formals (performance)$fuzzy), FALSE)
  expect_error (performance (NB, d [, -5], d [, 5], type = "roc", protocol = "crossvalidation",
                             nfolds = 3, seed = 0), NA)
  expect_error (performance (NB, d [, -5], d [, 5], type = "roc", fuzzy = FALSE,
                             protocol = "crossvalidation", nfolds = 3, seed = 0), NA)
  expect_error (performance (NB, d [, -5], d [, 5], type = "cost", protocol = "crossvalidation",
                             nfolds = 3, seed = 0), NA)
})

# --- Regression: BAGGING() and GRADIENTBOOSTING() rejected performance()'s arguments -------

test_that ("BAGGING() and GRADIENTBOOSTING() accept the standard method interface", {
  for (m in list (BAGGING, GRADIENTBOOSTING))
    expect_true (all (c ("tune", "methodparameters", "graph") %in% names (formals (m))))
  data (iris)
  expect_s3_class (BAGGING (iris [, -5], iris [, 5], LDA, nsamples = 3, tune = TRUE), "params")
})

# --- Regression: CDA() and LR() returned NULL on character labels --------------------------

test_that ("CDA() and LR() accept character class labels", {
  skip_if_not_installed ("MASS")
  data (iris)
  labels = as.character (iris [, 5])
  expect_s3_class (CDA (iris [, -5], labels), "cda")
  expect_s3_class (LR (iris [, -5], labels), "model")
})

test_that ("CDA() and LR() warn about empty classes instead of returning NULL", {
  skip_if_not_installed ("MASS")
  data (iris)
  d = iris [1:100, ]   # 'virginica' is a declared but unobserved level
  expect_warning (model <- CDA (d [, -5], d [, 5]), "no observation")
  expect_s3_class (model, "cda")
  expect_equal (model$nb.classes, 2)
  expect_warning (LR (d [, -5], d [, 5]), "no observation")
})

test_that ("CDA() fails clearly when there is only one class left", {
  data (iris)
  d = droplevels (iris [1:50, ])
  expect_error (CDA (d [, -5], d [, 5]), "at least 2 non-empty classes")
})

# --- Regression: STUMP() built its stump on a randomly drawn variable by default ------------

test_that ("STUMP() is deterministic by default and random on demand", {
  skip_if_not_installed ("rpart")
  data (iris)
  expect_equal (eval (formals (STUMP)$randomvar), FALSE)
  a = STUMP (iris [, -5], iris [, 5])
  b = STUMP (iris [, -5], iris [, 5])
  expect_equal (predict (a, iris [, -5]), predict (b, iris [, -5]))
  # The old behaviour is still reachable, and reproducible when seeded.
  expect_equal (predict (STUMP (iris [, -5], iris [, 5], randomvar = TRUE, seed = 0), iris [, -5]),
                predict (STUMP (iris [, -5], iris [, 5], randomvar = TRUE, seed = 0), iris [, -5]))
})

# --- Regression: sixteen S4 classes were declared but never instantiated --------------------

test_that ("the package's classes still answer is() and inherits() without setClass()", {
  data (iris)
  model = NB (iris [, -5], iris [, 5])
  expect_true (methods::is (model, "model"))
  expect_true (inherits (model, "model"))
  expect_true (methods::is (LDA (iris [, -5], iris [, 5]), "model"))
})

# =========================================================================================
# Ninth batch: multi-class measures, train/test evaluation, print methods
# =========================================================================================

# --- performance() on an explicit training set and test set --------------------------------

test_that ("performance() evaluates on a supplied test set, by default", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  model = NB (d$train.x, d$train.y)
  expected = evaluation.accuracy (predict (model, d$test.x), d$test.y)
  # No protocol named: the presence of a test set selects "holdout".
  expect_equal (unname (performance (NB, d$train.x, d$train.y, d$test.x, d$test.y)), expected)
  expect_equal (unname (performance (NB, d$train.x, d$train.y, d$test.x, d$test.y,
                                     protocol = "holdout")), expected)
})

test_that ("performance() refuses a protocol that would ignore the test set", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  for (p in c ("bootstrap", "crossvalidation", "loocv", "train"))
    expect_error (performance (NB, d$train.x, d$train.y, d$test.x, d$test.y, protocol = p),
                  "would ignore them", info = p)
  # Without a test set every protocol still works as before.
  expect_error (performance (NB, iris [, -5], iris [, 5], protocol = "bootstrap", nruns = 2,
                             seed = 0), NA)
})

# --- Multi-class precision, recall and everything derived from them -------------------------

test_that ("the pairwise measures work beyond two classes, and default to macro", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pred = predict (NB (d$train.x, d$train.y), d$test.x)
  for (f in list (evaluation.precision, evaluation.recall, evaluation.fmeasure,
                  evaluation.goodness, evaluation.jaccard, evaluation.fowlkesmallows))
  {
    expect_equal (f (pred, d$test.y), f (pred, d$test.y, average = "macro"))
    expect_length (f (pred, d$test.y), 1)
    expect_length (f (pred, d$test.y, average = "none"), 3)
    expect_named (f (pred, d$test.y, average = "none"), levels (d$test.y))
  }
})

test_that ("micro averaging coincides with accuracy", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pred = predict (NB (d$train.x, d$train.y), d$test.x)
  acc = evaluation.accuracy (pred, d$test.y)
  expect_equal (evaluation.precision (pred, d$test.y, average = "micro"), acc)
  expect_equal (evaluation.recall (pred, d$test.y, average = "micro"), acc)
  expect_equal (evaluation.fmeasure (pred, d$test.y, average = "micro"), acc)
})

test_that ("the per-class values match figures computed by hand", {
  # 3 classes, hand-built confusion matrix.
  gt   = factor (c ("a", "a", "a", "a", "b", "b", "b", "c", "c", "c"))
  pred = factor (c ("a", "a", "a", "b", "b", "b", "c", "c", "a", "b"))
  #        predicted a  b  c
  # true a           3  1  0
  # true b           0  2  1
  # true c           1  1  1
  expect_equal (unname (evaluation.precision (pred, gt, average = "none")),
                c (3 / 4, 2 / 4, 1 / 2))
  expect_equal (unname (evaluation.recall (pred, gt, average = "none")),
                c (3 / 4, 2 / 3, 1 / 3))
  expect_equal (evaluation.precision (pred, gt, average = "macro"),
                mean (c (3 / 4, 2 / 4, 1 / 2)))
  expect_equal (evaluation.recall (pred, gt, average = "weighted"),
                stats::weighted.mean (c (3 / 4, 2 / 3, 1 / 3), c (4, 3, 3)))
  expect_equal (evaluation.precision (pred, gt, average = "micro"), 6 / 10)
})

test_that ("two-class results are unchanged, and still driven by 'positive'", {
  gt   = factor (c ("+", "+", "+", "+", "-", "-", "-", "-", "-", "-"))
  pred = factor (c ("+", "+", "-", "-", "+", "-", "-", "-", "-", "-"))
  # TP = 2, FP = 1, FN = 2 for "+"
  expect_equal (evaluation.precision (pred, gt), 2 / 3)
  expect_equal (evaluation.recall (pred, gt), 2 / 4)
  expect_equal (evaluation.fmeasure (pred, gt), 2 * (2 / 3) * (1 / 2) / ((2 / 3) + (1 / 2)))
  # The other class as positive: TP = 5, FP = 2, FN = 1
  expect_equal (evaluation.precision (pred, gt, positive = "-"), 5 / 7)
  expect_equal (evaluation.recall (pred, gt, positive = "-"), 5 / 6)
  # Binary is the default on two classes, and macro is not the same thing.
  expect_equal (evaluation.precision (pred, gt),
                evaluation.precision (pred, gt, average = "binary"))
  expect_false (isTRUE (all.equal (evaluation.precision (pred, gt),
                                   evaluation.precision (pred, gt, average = "macro"))))
})

test_that ("the measures explain themselves when they cannot answer", {
  data (iris)
  pred = iris [, 5]
  expect_error (evaluation.precision (pred, iris [, 5], average = "binary"), "two classes")
  expect_error (evaluation.precision (pred, iris [, 5], average = "banana"), "should be one of")
  gt = factor (c ("+", "+", "-", "-"))
  expect_error (evaluation.precision (gt, gt, positive = "?"), "not one of the class")
  # average = "none" cannot go through evaluation(), which returns one number per criterion.
  expect_error (evaluation (pred, iris [, 5], eval = "precision", average = "none"),
                "returned 3 values")
})

test_that ("a class nobody predicts scores 0 rather than NaN", {
  gt   = factor (c ("a", "a", "b", "b", "c", "c"))
  pred = factor (c ("a", "a", "b", "b", "b", "b"), levels = c ("a", "b", "c"))
  p = evaluation.precision (pred, gt, average = "none")
  expect_equal (unname (p ["c"]), 0)
  expect_false (any (is.nan (p)))
  expect_false (any (is.nan (evaluation.fmeasure (pred, gt, average = "none"))))
})

test_that ("evaluation() and performance() carry the multi-class criteria through", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pred = predict (NB (d$train.x, d$train.y), d$test.x)
  res = evaluation (pred, d$test.y, eval = c ("accuracy", "precision", "recall", "fmeasure"))
  expect_named (res, c ("accuracy", "precision", "recall", "fmeasure"))
  expect_true (all (res > 0.5))
  perf = performance (c (NB, LDA), d$train.x, d$train.y, d$test.x, d$test.y,
                      eval = c ("precision", "recall"))
  expect_equal (dim (perf), c (2, 2))
})

# --- print methods --------------------------------------------------------------------------

test_that ("the print methods say what the object is instead of dumping it", {
  data (iris)
  shown = function (x) paste (utils::capture.output (print (x)), collapse = "\n")
  expect_match (shown (NB (iris [, -5], iris [, 5])), "Classification/regression model \\(NB\\)")
  expect_match (shown (splitdata (iris, 5, seed = 0)), "Training/test split")
  expect_match (shown (splitdata (iris, 5, seed = 0)), "105 observations")
  expect_match (shown (KNN (iris [, -5], iris [, 5])), "K-nearest-neighbours")
  expect_match (shown (CDA (iris [, -5], iris [, 5])), "Canonical discriminant analysis")
  expect_match (shown (BAGGING (iris [, -5], iris [, 5], LDA, nsamples = 3)), "Ensemble model")
  expect_match (shown (SVM (iris [, -5], iris [, 5], tune = TRUE)), "gamma")
  expect_match (shown (NB (iris [, -5], iris [, 5], tune = TRUE)), "none")
  sel = selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "cfs")
  expect_match (shown (sel), "Feature selection")
  # The names of the selected features, not only their positions.
  expect_match (shown (sel), "Petal")
  expect_equal (sel$features, colnames (iris [, -5]) [sel$selection])
  # Short: the whole point is not to flood the console.
  expect_lt (length (utils::capture.output (print (NB (iris [, -5], iris [, 5])))), 6)
})

test_that ("summary() of a model still gives the detail of the wrapped model", {
  data (iris)
  out = utils::capture.output (summary (NB (iris [, -5], iris [, 5])))
  expect_match (paste (out, collapse = "\n"), "Classification/regression model")
  expect_gt (length (out), 5)
})

# =========================================================================================
# Tenth batch: penalized logistic regression
# =========================================================================================

test_that ("LR() keeps its unpenalized behaviour by default", {
  data (iris)
  a = LR (iris [, -5], iris [, 5])
  expect_equal (a$method, "LR")
  expect_equal (eval (formals (LR)$reg) [1], "none")
})

test_that ("LR() fits the three penalties, on two classes and on more", {
  skip_if_not_installed ("glmnet")
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  for (r in c ("ridge", "lasso", "elastic"))
  {
    model = suppressMessages (LR (d$train.x, d$train.y, reg = r))
    expect_equal (model$method, "penalizedlr", info = r)
    expect_true (model$lambda > 0)
    pred = predict (model, d$test.x)
    expect_s3_class (pred, "factor")
    expect_equal (levels (pred), levels (d$train.y))
    expect_gt (evaluation.accuracy (pred, d$test.y), 0.8)
  }
  # Two classes go through glmnet's binomial family, which returns a different shape.
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  s = splitdata (two, 5, seed = 0)
  model = suppressMessages (LR (s$train.x, s$train.y, reg = "ridge"))
  expect_equal (levels (predict (model, s$test.x)), levels (s$train.y))
})

test_that ("the probabilities of a penalized LR agree with its hard predictions", {
  # The bug of the first batch (predictmodel.lr, columns swapped) in a new implementation:
  # argmax of the probabilities must be the predicted class, on 2 classes and on 3.
  skip_if_not_installed ("glmnet")
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  for (dataset in list (d, { two = iris
                             levels (two [, 5]) = c ("+", "+", "-")
                             splitdata (two, 5, seed = 0) }))
  {
    model = suppressMessages (LR (dataset$train.x, dataset$train.y, reg = "lasso"))
    hard = predict (model, dataset$test.x)
    fuzzy = predict (model, dataset$test.x, fuzzy = TRUE)
    expect_equal (colnames (fuzzy), levels (dataset$train.y))
    expect_equal (colnames (fuzzy) [apply (fuzzy, 1, which.max)], as.character (hard))
    expect_equal (unname (rowSums (fuzzy)), rep (1, nrow (fuzzy)), tolerance = 1e-8)
  }
})

test_that ("a penalized LR is usable through performance(), including ROC curves", {
  skip_if_not_installed ("glmnet")
  skip_if_not_installed ("ROCR")
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  d = splitdata (two, 5, seed = 0)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (suppressMessages (
    performance (LR, d$train.x, d$train.y, d$test.x, d$test.y, reg = "ridge")), NA)
  expect_error (suppressMessages (
    performance (LR, d$train.x, d$train.y, d$test.x, d$test.y, type = "roc", reg = "ridge")), NA)
})

# =========================================================================================
# Eleventh batch: boosting
# =========================================================================================

separable <-
  function ()
  {
    set.seed (0)
    list (x = data.frame (V1 = c (stats::rnorm (50, 0), stats::rnorm (50, 20)),
                          V2 = c (stats::rnorm (50, 0), stats::rnorm (50, 20))),
          y = factor (rep (c ("a", "b"), each = 50)))
  }

test_that ("ADABOOST() keeps a base learner that is perfect on the training set", {
  d = separable ()
  # LDA classifies these two clouds without a mistake, so epsilon is 0 (or 1e-16). That model
  # used to fall through the 'if (epsilon > 0)' branch: it was discarded, the weights were not
  # updated, and the loop refitted and discarded the same model 'nsamples' times. ADABOOST
  # returned an ensemble of *zero* models and predict() then failed with
  # "dim(X) must have a positive length".
  for (fuzzy in c (FALSE, TRUE))
  {
    model = suppressMessages (ADABOOST (d$x, d$y, LDA, nsamples = 100, fuzzy = fuzzy))
    expect_gte (length (model$models), 1)
    pred = predict (model, d$x)
    expect_length (pred, nrow (d$x))
    expect_equal (as.character (pred), as.character (d$y))
  }
})

test_that ("ADABOOST() records how many models it actually kept", {
  data (iris)
  # It used to announce it with a message on every call; the number is on the object, and
  # print.boosting() is where it is said. See the test at the end of this file.
  expect_message (model <- ADABOOST (iris [, -5], iris [, 5], LDA, nsamples = 100), NA)
  expect_lt (length (model$models), 100)
  expect_equal (model$nsamples, 100)
  # ... and an ensemble that kept them all is silent too.
  expect_silent (BAGGING (iris [, -5], iris [, 5], LDA, nsamples = 5))
})

test_that ("an ensemble that could keep nothing fails with an explanation", {
  data (iris)
  # A 'learning method' whose predictions are pure noise cannot beat chance.
  useless = function (train, labels, ...)
  {
    res = list (lev = levels (factor (labels)), method = "useless")
    class (res) = "useless"
    return (res)
  }
  registerS3method ("predict", "useless",
                    function (object, test, fuzzy = FALSE, ...)
                      factor (rep (object$lev [1], nrow (as.data.frame (test))),
                              levels = object$lev),
                    envir = environment ())
  expect_error (suppressMessages (ADABOOST (iris [, -5], iris [, 5], useless, nsamples = 5)),
                "not a single model")
})

test_that ("predict() on an ensemble weights the models in fuzzy mode too", {
  data (iris)
  model = suppressMessages (ADABOOST (iris [, -5], iris [, 5], STUMP, nsamples = 20))
  weights = sapply (model$models, function (m) m$boostweight)
  skip_if (length (unique (round (weights, 8))) < 2)
  fuzzy = predict (model, iris [, -5], fuzzy = TRUE)
  # The plain (unweighted) mean the previous version computed.
  unweighted = Reduce ("+", lapply (model$models,
                                    function (m) as.matrix (predict (m, iris [, -5],
                                                                     fuzzy = TRUE)))) /
               length (model$models)
  expect_false (isTRUE (all.equal (unname (fuzzy), unname (unweighted))))
  expect_equal (unname (rowSums (fuzzy)), rep (1, nrow (iris)), tolerance = 1e-8)
  expect_equal (colnames (fuzzy), levels (iris [, 5]))
})

test_that ("the ensemble vote is unchanged, and fast", {
  data (iris)
  model = BAGGING (iris [, -5], iris [, 5], LDA, nsamples = 20)
  pred = predict (model, iris [, -5])
  # All the weights of a bagging are 1, so the vote is the plain majority; computed here the
  # slow way the previous version did.
  raw = sapply (model$models, function (m) as.character (predict (m, iris [, -5])))
  majority = apply (raw, 1, function (v) names (which.max (table (v))))
  expect_equal (as.character (pred), majority)
  expect_s3_class (pred, "factor")
  expect_equal (levels (pred), levels (iris [, 5]))
})

# =========================================================================================
# Second audit, B1 / B4
# =========================================================================================

test_that ("CDA reports the share of the trace, not the share of the squares", {
  data (iris)
  model = CDA (iris [, -5], iris [, 5])
  eig = model$eig
  # 'percentage of variance' is the share of the trace of the eigenvalues shown next to it.
  # It used to be computed on their squares, which is the share of nothing.
  expect_equal (unname (eig [, 2]), unname (100 * eig [, 1] / sum (eig [, 1])))
  expect_equal (unname (eig [, 3]), unname (cumsum (eig [, 2])))
  expect_equal (unname (eig [nrow (eig), 3]), 100)
})

test_that ("CDA's discriminant power is the proportion of trace of MASS::lda", {
  skip_if_not_installed ("MASS")
  data (iris)
  # CDA diagonalises V^-1 B, MASS::lda works on W^-1 B; the two describe the same axes and
  # lambda_W = lambda_V / (1 - lambda_V), so the shares must agree exactly.
  model = CDA (iris [, -5], iris [, 5])
  reference = MASS::lda (iris [, -5], iris [, 5])
  expect_equal (unname (model$eig [, 4]),
                unname (100 * reference$svd^2 / sum (reference$svd^2)))
})

test_that ("CDA's eigenvalue table survives a single axis and a perfect separation", {
  data (iris)
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  eig = CDA (two [, -5], two [, 5])$eig
  expect_null (dim (eig))
  expect_named (eig, c ("eigenvalue", "percentage of variance",
                        "cumulative percentage of variance", "discriminant power"))
  expect_equal (unname (eig [2:4]), c (100, 100, 100))
  # A perfectly separating axis has lambda = 1, hence infinite discriminant power: the whole
  # of it goes to that axis rather than to NaN.
  x = data.frame (a = c (seq (0, 1, length.out = 20), seq (50, 51, length.out = 20)),
                  b = c (1:20, 20:1))
  y = factor (rep (c ("A", "B"), each = 20))
  expect_true (all (is.finite (CDA (x, y)$eig)))
  expect_equal (fdm2id:::cda.power (c (1, 0.5)), c (1, 0))
  expect_equal (sum (fdm2id:::cda.power (c (0.9, 0.2))), 1)
})

test_that ("kappa puts predictions and ground truth on a common set of labels", {
  # as.numeric() on two factors with different levels compares codes that do not stand for
  # the same classes -- which is what happens whenever a model never predicts some class.
  gt = factor (c (rep ("a", 5), rep ("b", 5), rep ("c", 5)))
  short = factor (c (rep ("a", 5), rep ("c", 5), rep ("c", 5)), levels = c ("a", "c"))
  full = factor (as.character (short), levels = c ("a", "b", "c"))
  expect_equal (evaluation.kappa (short, gt), evaluation.kappa (full, gt))
})

test_that ("kappa is Cohen's unweighted kappa, so class labels stay nominal", {
  # weight = "equal" is a *linearly weighted* kappa: it reads the class codes as an ordinal
  # scale, and the same number of errors then scores differently depending on which classes
  # were confused.
  gt = factor (rep (c ("a", "b", "c"), each = 20))
  far = gt; far [1:10] = "c"
  near = gt; near [1:10] = "b"
  expect_equal (evaluation.kappa (far, gt), evaluation.kappa (near, gt))
  # The value itself, from the definition.
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pred = predict (NB (d$train.x, d$train.y), d$test.x)
  t = table (d$test.y, pred)
  n = sum (t)
  po = sum (diag (t)) / n
  pe = sum (rowSums (t) * colSums (t)) / (n * n)
  expect_equal (evaluation.kappa (pred, d$test.y), (po - pe) / (1 - pe))
  # One class, always right: chance agreement is already perfect and kappa is 0 / 0.
  expect_true (is.na (evaluation.kappa (factor (rep ("a", 10)), factor (rep ("a", 10)))))
})

# --- Fifth audit: performance() always names its rows --------------------------------------

test_that ("performance() names its methods even when they come from a variable", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  expect_equal (rownames (performance (c (NB, LDA), d$train.x, d$train.y,
                                       d$test.x, d$test.y)), c ("NB", "LDA"))
  # A list held in a variable carries no names in the call, and the rows used to come back
  # unnamed rather than numbered.
  methods = c (NB, LDA)
  expect_equal (rownames (performance (methods, d$train.x, d$train.y, d$test.x, d$test.y)),
                c ("Method 1", "Method 2"))
  expect_equal (rownames (performance (methods, d$train.x, d$train.y, d$test.x, d$test.y,
                                       names = c ("a", "b"))), c ("a", "b"))
})

# =========================================================================================
# Eighth audit: plot.cda
# =========================================================================================

test_that ("plot.cda draws on the whole device and leaves it as it found it", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  d = splitdata (iris, 5, seed = 0)
  model = CDA (d$train.x, d$train.y)
  # The legend used to sit in a band of its own, obtained with layout(): an eighth of the
  # device whatever its size, which truncated the entries on a wide screen.
  before = graphics::par (c ("mfrow", "mar", "mfg"))
  expect_error (plot (model), NA)
  expect_equal (graphics::par (c ("mfrow", "mar", "mfg")), before)
  expect_error (plot (model, d$test.x), NA)
  expect_equal (graphics::par (c ("mfrow", "mar", "mfg")), before)
  expect_error (plot (model, d$test.x, legendpos = "bottomright"), NA)
  expect_error (plot (model, d$test.x, axes = c (2, 1)), NA)
})

test_that ("plot.cda works on a single canonical axis", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  d = splitdata (two, 5, seed = 0)
  model = CDA (d$train.x, d$train.y)
  expect_equal (model$dim, 1)
  expect_error (plot (model), NA)
  # 'n' was the number of training observations, computed before the test set was appended,
  # so cbind (1:n, X) recycled and the test points were drawn at the wrong index.
  expect_warning (plot (model, d$test.x), NA)
  expect_error (plot (model, d$test.x), NA)
})

test_that ("a canonical axis is labelled with the variance it carries", {
  data (iris)
  model = CDA (iris [, -5], iris [, 5])
  expect_equal (fdm2id:::cda.axislabel (model, 1),
                paste0 ("Can. 1 (", round (model$eig [1, 2], 2), " %)"))
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  flat = CDA (two [, -5], two [, 5])
  # With a single axis the table is a named vector, not a matrix.
  expect_equal (fdm2id:::cda.axislabel (flat, 1), "Can. 1 (100 %)")
})

# --- Tenth audit: an ensemble says its size instead of announcing it ------------------------

test_that ("ADABOOST stops early without a message, and says so when printed", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  # The message fired on every call; the information belongs to whoever looks at the model.
  expect_message (model <- ADABOOST (d$train.x, d$train.y, STUMP, seed = 0), NA)
  expect_lt (length (model$models), 100)
  expect_equal (model$nsamples, 100)
  expect_output (print (model), "of the 100 asked for")
  # An ensemble that keeps everything says nothing of the sort.
  full = BAGGING (d$train.x, d$train.y, NB, nsamples = 5, seed = 0)
  expect_length (full$models, 5)
  expect_output (print (full), "models")
})

# --- Regression: predict.model() ignored fuzzy for RANDOMFOREST -------------------------
# randomForest() needs type = "prob" to return scores; without it fuzzy = TRUE gave back hard
# labels, and performance (type = "roc") then had no column for the positive class.

test_that ("predict (RANDOMFOREST, fuzzy = TRUE) returns one column per class", {
  skip_if_not_installed ("randomForest")
  data (iris)
  m = RANDOMFOREST (iris [, -5], iris [, 5], seed = 0)
  p = predict (m, iris [, -5], fuzzy = TRUE)
  expect_true (is.matrix (p))
  expect_equal (colnames (p), levels (iris [, 5]))
  expect_equal (nrow (p), nrow (iris))
  expect_equal (unname (rowSums (p)), rep (1, nrow (iris)))
})

test_that ("performance (type = 'roc') works with RANDOMFOREST", {
  skip_if_not_installed ("randomForest")
  skip_if_not_installed ("ROCR")
  data (ionosphere)
  expect_error (
    performance (RANDOMFOREST, ionosphere [, -34], ionosphere [, 34], type = "roc",
                 protocol = "holdout", seed = 0),
    NA
  )
})

# --- Regression: confusion() rejected the ... that performance() forwards ----------------
# performance() hands its own ... both to the learning method and to the evaluation function,
# so an ensemble's learningmethod reached confusion() as an unused argument.

test_that ("performance (type = 'confusion') accepts a learningmethod", {
  skip_if_not_installed ("rpart")
  data (iris)
  expect_error (
    performance (BAGGING, iris [, -5], iris [, 5], type = "confusion", protocol = "holdout",
                 seed = 0, learningmethod = CART, graph = FALSE),
    NA
  )
})

# --- Regression: cbind() of factors reached the curves as integer level codes -------------
# cbind (predict (m1, x), predict (m2, x)) drops the factor class. Read as scores, the codes
# rank the observations by level number: the curve came out upside down whenever 'positive'
# was not the last level, and type = "hard" compared "1"/"2" against the labels and scored
# everything 0.

test_that ("curve.scores() maps a cbind() of factors back to the class labels", {
  data (iris)
  d = iris
  levels (d [, 5]) = c ("+", "+", "-")
  hard = factor (rep (c ("+", "-"), length.out = nrow (d)), levels = c ("+", "-"))
  m = cbind (hard, hard)
  for (type in c ("auto", "hard"))
  {
    s = fdm2id:::curve.scores (m, d [, 5], "+", type)
    expect_equal (s [[1]], as.numeric (hard == "+"), info = type)
    expect_equal (s [[2]], as.numeric (hard == "+"), info = type)
  }
})

test_that ("roc.curves() gives the same curve whichever level is the positive class", {
  skip_if_not_installed ("ROCR")
  skip_if_not_installed ("e1071")
  data (iris)
  d = iris
  levels (d [, 5]) = c ("+", "+", "-")
  model = NB (d [, -5], d [, 5])
  p = cbind (predict (model, d [, -5]), predict (model, d [, -5]))
  auc = function (positive)
  {
    s = fdm2id:::curve.scores (p, d [, 5], positive, "auto") [[1]]
    lab = as.numeric (d [, 5] == positive)
    unlist (ROCR::performance (ROCR::prediction (s, lab), "auc")@y.values)
  }
  # The same predictions, read from either side, are equally good: both curves sit above the
  # diagonal. Reading the level codes as scores put one of the two below it.
  expect_gt (auc ("+"), 0.5)
  expect_gt (auc ("-"), 0.5)
})
