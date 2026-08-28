# =========================================================================================
# The uniform contract the package promises
#
# fdm2id's whole point is that every method is called the same way. These tests check that
# claim mechanically rather than by reading the signatures, so that a new method cannot join
# the package with a different convention without the suite noticing.
# =========================================================================================

supervised <-
  function ()
    c ("APRIORI", "CART", "CDA", "KNN", "LDA", "LR", "MLP", "NB", "QDA", "RANDOMFOREST",
       "STUMP", "SVM", "SVMl", "SVMr", "GRADIENTBOOSTING", "ADABOOST", "BAGGING",
       "FEATURESELECTION",
       "GBREG", "KERREG", "LINREG", "MLPREG", "POLYREG", "SVR", "SVRl", "SVRr")

clusterings <-
  function ()
    c ("DBSCAN", "EM", "HCA", "KMEANS", "MEANSHIFT", "PAM", "SOM", "SPECTRAL")

test_that ("every learning method takes the same first two arguments", {
  for (name in supervised ())
  {
    args = names (formals (get (name)))
    expect_true (identical (args [1:2], c ("train", "labels")) ||
                 identical (args [1:2], c ("x", "y")),
                 info = paste (name, ":", paste (args [1:2], collapse = ", ")))
  }
})

test_that ("every learning method ends on the same block of arguments", {
  # tune, methodparameters, graph and seed are what performance() passes to whatever method it
  # is given, so all of them must accept all four, in the same order, and then '...'. A method
  # that has nothing to draw at random still takes 'seed', so that a script can be moved from
  # one method to another without rewriting the call.
  for (name in supervised ())
  {
    args = names (formals (get (name)))
    expect_identical (tail (args, 5),
                      c ("tune", "methodparameters", "graph", "seed", "..."),
                      info = paste (name, ":", paste (args, collapse = ", ")))
  }
})

generators <-
  function ()
    c ("data.diag", "data.gauss", "data.parabol", "data.target1", "data.target2",
       "data.twomoons", "data.xor")

test_that ("no learning method draws unless it is asked to", {
  # A function that plots by default cannot be called in a loop, a report or a vignette
  # without a device piling up. The seven data generators were the last ones that did.
  for (name in c (supervised (), clusterings (), generators ()))
  {
    f = formals (get (name))
    if (!("graph" %in% names (f)))
      next
    expect_identical (f$graph, FALSE, info = name)
  }
})

test_that ("every clustering takes its number of clusters second", {
  # HCA () read its second argument as the linkage, so HCA (iris [, -5], 3) silently asked
  # for a method named "3" instead of three clusters.
  for (name in c ("EM", "HCA", "KMEANS", "PAM", "SPECTRAL"))
  {
    args = names (formals (get (name)))
    expect_identical (args [1:2], c ("d", "k"),
                      info = paste (name, ":", paste (args [1:2], collapse = ", ")))
  }
  # The density-based ones have no number of clusters to take: it is what they find.
  for (name in c ("DBSCAN", "MEANSHIFT"))
    expect_false ("k" %in% names (formals (get (name))) [1:2], info = name)
})

test_that ("every learning method answers tune = TRUE with a params object", {
  data (iris)
  data (trees)
  regression = c ("GBREG", "KERREG", "LINREG", "MLPREG", "POLYREG", "SVR", "SVRl", "SVRr")
  for (name in setdiff (supervised (), c ("ADABOOST", "BAGGING", "FEATURESELECTION")))
  {
    x = if (name %in% regression) trees [, -3] else iris [, -5]
    y = if (name %in% regression) trees [, 3] else iris [, 5]
    res = suppressWarnings (suppressMessages (do.call (name, list (x, y, tune = TRUE))))
    expect_s3_class (res, "params")
  }
})

test_that ("a method with nothing to tune answers with an empty params object", {
  # Those methods expose no hyperparameter, so tune = TRUE is answered without fitting
  # anything -- and the empty object they return must leave every default alone when
  # performance() hands it back (see the last test of this file).
  data (iris)
  for (name in c ("CART", "CDA", "LDA", "NB", "QDA", "RANDOMFOREST", "STUMP",
                  "GRADIENTBOOSTING", "APRIORI"))
  {
    res = suppressWarnings (suppressMessages (
      do.call (name, list (iris [, -5], iris [, 5], tune = TRUE))))
    expect_s3_class (res, "params")
    expect_length (res, 0)
  }
})

test_that ("a method that tunes answers with the parameters it retained, and takes them back", {
  # This is what tune = TRUE is for: performance() asks once, then hands the answer back on
  # every split instead of re-running the search. MLPREG() returned an empty object, so its
  # network was re-tuned on every split and it was the one method of the pair MLP/MLPREG that
  # could not be pre-tuned.
  data (iris)
  data (trees)
  cases = list (list (MLP, iris [, -5], iris [, 5], c ("decay", "hidden")),
                list (MLPREG, trees [, -3], trees [, 3], c ("decay", "hidden")),
                list (SVM, iris [, -5], iris [, 5], c ("gamma", "cost")),
                list (SVR, trees [, -3], trees [, 3], c ("epsilon", "gamma", "cost")))
  for (case in cases)
  {
    method = case [[1]]
    params = suppressWarnings (suppressMessages (
      method (case [[2]], case [[3]], nfolds = 3, tune = TRUE, seed = 0)))
    expect_s3_class (params, "params")
    expect_named (params, case [[4]])
    expect_true (all (sapply (params, length) == 1))
    # And the model fitted from them is a model, not another round of tuning.
    expect_s3_class (suppressWarnings (suppressMessages (
      method (case [[2]], case [[3]], methodparameters = params, seed = 0))), "model")
  }
})

test_that ("FEATURESELECTION answers tune = TRUE like the method it wraps", {
  # It used to run the whole selection and return a "selection" object whose $model was a
  # params one -- neither what performance() asks for nor anything it can hand back.
  data (iris)
  expect_s3_class (FEATURESELECTION (iris [, -5], iris [, 5], uninb = 2, mainmethod = LDA,
                                     tune = TRUE), "params")
  tuned = suppressWarnings (suppressMessages (
    FEATURESELECTION (iris [, -5], iris [, 5], uninb = 2, mainmethod = SVM,
                      gamma = 2^(0:1), cost = 2^(0:1), nfolds = 3, tune = TRUE, seed = 0)))
  expect_s3_class (tuned, "params")
  expect_named (tuned, c ("gamma", "cost"))
})

test_that ("every learning method can be handed to performance()", {
  skip_on_cran ()
  data (iris)
  data (trees)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  classification = c ("CART", "CDA", "KNN", "LDA", "LR", "MLP", "NB", "QDA", "STUMP",
                      "SVM", "SVMl", "SVMr")
  d = splitdata (iris, 5, seed = 0)
  for (name in classification)
  {
    res = suppressWarnings (suppressMessages (
      performance (get (name), d$train.x, d$train.y, d$test.x, d$test.y)))
    expect_length (res, 1)
    expect_true (is.finite (res), info = name)
  }
  regression = c ("GBREG", "KERREG", "LINREG", "MLPREG", "POLYREG", "SVR", "SVRl", "SVRr")
  r = splitdata (trees, 3, seed = 0)
  for (name in regression)
  {
    res = suppressWarnings (suppressMessages (
      performance (get (name), r$train.x, r$train.y, r$test.x, r$test.y)))
    expect_length (res, 1)
    expect_true (is.finite (res), info = name)
  }
})

test_that ("performance() accepts an explicit 'graph', rather than clashing with its own", {
  data (trees)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  r = splitdata (trees, 3, seed = 0)
  # LINREG draws by default; silencing it inside performance() used to stop on
  # 'formal argument "graph" matched by multiple actual arguments'.
  expect_error (performance (LINREG, r$train.x, r$train.y, r$test.x, r$test.y, graph = FALSE),
                NA)
  expect_equal (performance (LINREG, r$train.x, r$train.y, r$test.x, r$test.y, graph = FALSE),
                performance (LINREG, r$train.x, r$train.y, r$test.x, r$test.y))
})

test_that ("a method that draws at random is reproducible when seeded", {
  data (iris)
  data (trees)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  twice = function (name, x, y)
  {
    a = predict (do.call (name, list (x, y, seed = 0)), x)
    b = predict (do.call (name, list (x, y, seed = 0)), x)
    expect_equal (as.character (a), as.character (b), info = name)
  }
  for (name in c ("KNN", "SVM", "MLP", "CART", "RANDOMFOREST", "LR"))
    twice (name, iris [, -5], iris [, 5])
  for (name in c ("MLPREG", "SVR", "LINREG"))
    suppressWarnings (twice (name, trees [, -3], trees [, 3]))
})

test_that ("every clustering returns its assignments in the same place", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  x = iris [, -5]
  built = list (DBSCAN = suppressMessages (DBSCAN (x)), EM = EM (x, 3), HCA = HCA (x, 3),
                KMEANS = KMEANS (x, k = 3, seed = 0),
                MEANSHIFT = MEANSHIFT (x), PAM = PAM (x, 3), SOM = SOM (x, 4, 4),
                SPECTRAL = SPECTRAL (x, 3))
  expect_named (built, clusterings ())
  for (name in names (built))
  {
    res = built [[name]]
    expect_false (is.null (res$cluster), info = name)
    expect_length (res$cluster, nrow (x))
    expect_true (is.numeric (res$cluster) || is.integer (res$cluster), info = name)
  }
})

test_that ("every clustering can place a new observation", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  d = splitdata (iris, 5, seed = 0)
  built = list (DBSCAN = suppressMessages (DBSCAN (d$train.x)), EM = EM (d$train.x, 3),
                HCA = HCA (d$train.x, 3),
                KMEANS = KMEANS (d$train.x, k = 3, seed = 0), MEANSHIFT = MEANSHIFT (d$train.x),
                PAM = PAM (d$train.x, 3), SOM = SOM (d$train.x, 4, 4),
                SPECTRAL = SPECTRAL (d$train.x, 3))
  for (name in names (built))
  {
    res = predict (built [[name]], d$test.x)
    expect_length (res, nrow (d$test.x))
    expect_false (anyNA (res), info = name)
  }
})

test_that ("every model of the package answers predict() the same way", {
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  for (name in c ("NB", "LDA", "QDA", "CDA", "CART", "KNN", "LR", "SVM"))
  {
    model = do.call (name, list (d$train.x, d$train.y))
    hard = predict (model, d$test.x)
    expect_length (hard, nrow (d$test.x))
    expect_equal (levels (factor (hard)), levels (droplevels (factor (hard))), info = name)
    fuzzy = predict (model, d$test.x, fuzzy = TRUE)
    expect_equal (nrow (fuzzy), nrow (d$test.x), info = name)
    expect_equal (colnames (fuzzy), levels (d$train.y), info = name)
    # The hard prediction is the most likely class -- the first batch of the audit found two
    # methods where it was not.
    expect_equal (colnames (fuzzy) [apply (fuzzy, 1, which.max)], as.character (hard),
                  info = name)
  }
})

test_that ("a method refuses a target it cannot handle, saying which one to use", {
  data (iris)
  data (trees)
  expect_error (ADABOOST (trees [, -3], trees [, 3], LINREG, nsamples = 3), "BAGGING")
  expect_error (GBREG (iris [, -5], iris [, 5]), "GRADIENTBOOSTING")
})

test_that ("an empty 'params' object leaves the defaults alone", {
  # performance() asks each method for its hyperparameters (tune = TRUE) and hands the answer
  # back when fitting. A method with nothing to tune answers with an *empty* params object,
  # which must not overwrite anything: MLPREG() read size and decay out of it unconditionally,
  # got NULL, and stopped on "NA/NaN argument" inside nnet().
  data (iris)
  data (trees)
  empty = NB (iris [, -5], iris [, 5], tune = TRUE)
  expect_s3_class (empty, "params")
  expect_length (empty, 0)
  expect_s3_class (MLPREG (trees [, -3], trees [, 3], methodparameters = empty), "model")
  expect_s3_class (MLP (iris [, -5], iris [, 5], methodparameters = empty), "model")
  expect_s3_class (SVM (iris [, -5], iris [, 5], methodparameters = empty), "model")
  expect_s3_class (SVR (trees [, -3], trees [, 3], methodparameters = empty), "model")
})

# =========================================================================================
# Third audit, C1 -- the resampling scheme behind hyperparameter search
# =========================================================================================

test_that ("every method that searches a grid takes the same 'nfolds' argument", {
  # LR() already chose its penalty by an nfolds-fold cross-validation; the methods that search
  # a grid through e1071 used a fixed twenty-sample bootstrap nobody could reach.
  for (name in c ("KNN", "MLP", "SVM", "SVMl", "SVMr", "MLPREG", "SVR", "SVRl", "SVRr", "LR"))
  {
    args = names (formals (get (name)))
    expect_true ("nfolds" %in% args, info = name)
    # ... and it still ends on the common block.
    expect_identical (tail (args, 5), c ("tune", "methodparameters", "graph", "seed", "..."),
                      info = name)
  }
})

test_that ("'nfolds' reaches the tuning of the wrappers as well as of SVM/SVR themselves", {
  data (iris)
  data (trees)
  # A number of folds larger than the sample is what e1071 refuses, so it is a cheap way of
  # checking that the argument is actually forwarded rather than silently dropped.
  expect_error (fdm2id:::tune.scheme (1), "at least|>= 2")
  expect_error (fdm2id:::tune.scheme (NA), ">= 2")
  expect_s3_class (fdm2id:::tune.scheme (5), "tune.control")
  expect_equal (fdm2id:::tune.scheme (5)$cross, 5)
  expect_equal (fdm2id:::tune.scheme (5)$sampling, "cross")
  # A small grid and few folds, so the test stays fast.
  expect_s3_class (SVMl (iris [, -5], iris [, 5], cost = 2^(0:1), nfolds = 3, seed = 0), "model")
  expect_s3_class (SVRr (trees [, -3], trees [, 3], gamma = 2^(0:1), cost = 1, epsilon = .1,
                         nfolds = 3, seed = 0), "model")
})

test_that ("a tuned SVM still carries the probability model its predictions need", {
  # The grid is searched without probability = TRUE (e1071 fits that by an internal five-fold
  # cross-validation, six fits where one is needed, all thrown away), so the retained model
  # must be refitted with it.
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  model = SVM (d$train.x, d$train.y, gamma = 2^(0:1), cost = 2^(0:1), nfolds = 3, seed = 0)
  fuzzy = predict (model, d$test.x, fuzzy = TRUE)
  expect_equal (nrow (fuzzy), nrow (d$test.x))
  expect_setequal (colnames (fuzzy), levels (d$train.y))
  expect_true (all (abs (rowSums (fuzzy) - 1) < 1e-8))
  # The hard prediction is the most likely class.
  expect_equal (colnames (fuzzy) [apply (fuzzy, 1, which.max)],
                as.character (predict (model, d$test.x)))
})
