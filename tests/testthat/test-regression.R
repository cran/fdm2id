# Tests for regression.R

test_that ("LINREG smoke test: default linear fit", {
  data (trees)
  model = LINREG (trees [, -3], trees [, 3])
  expect_s3_class (model, "model")
  pred = predict (model, trees [, -3])
  expect_equal (length (pred), nrow (trees))
})

test_that ("LINREG's quali parameter (renamed from the former formula=) handles qualitative variables", {
  data (ToothGrowth)
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "intercept"), NA)
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "slope"), NA)
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "both"), NA)
})

test_that ("POLYREG / KERREG smoke tests", {
  skip_if_not_installed ("mda")
  skip_if_not_installed ("ibr")
  data (trees)
  poly = POLYREG (trees [, -3], trees [, 3])
  expect_equal (length (predict (poly, trees [, -3])), nrow (trees))
  ker = KERREG (trees [, -3], trees [, 3])
  expect_equal (length (predict (ker, trees [, -3])), nrow (trees))
})

test_that ("SVR smoke test", {
  skip_if_not_installed ("e1071")
  data (trees)
  model = SVR (trees [, -3], trees [, 3], kernel = "linear", cost = 1)
  expect_equal (length (predict (model, trees [, -3])), nrow (trees))
})

test_that ("evaluation.r2 / evaluation.msep return sensible numeric values", {
  data (trees)
  d = splitdata (trees, 3, seed = 0)
  model = LINREG (d$train.x, d$train.y)
  pred = predict (model, d$test.x)
  expect_true (is.numeric (evaluation.r2 (pred, d$test.y)))
  expect_true (is.numeric (evaluation.msep (pred, d$test.y)))
  expect_gte (evaluation.msep (pred, d$test.y), 0)
})

# --- Regression: LINREG() hard-coded the length of the default lambda grid ----------------
# lambda [101:1] reversed the *default* 101-value grid; any other grid gave NAs, and the
# model was then fitted on lambda = NA.

test_that ("LINREG() accepts a lambda grid of any length", {
  skip_if_not_installed ("glmnet")
  data (mtcars)
  grid = 10^seq (-3, 3, length.out = 25)
  for (r in c ("ridge", "lasso", "elastic"))
  {
    model = LINREG (mtcars [, -1], mtcars [, 1], reg = r, lambda = grid, graph = FALSE, nrep = 2)
    expect_false (is.na (model$model$lambda), info = r)
    # glmnet gives back the value it was handed, up to a floating-point round trip.
    expect_lt (min (abs (grid - model$model$lambda)), 1e-8)
  }
  expect_equal (length (predict (LINREG (mtcars [, -1], mtcars [, 1], reg = "ridge",
                                         lambda = grid, graph = FALSE, nrep = 2),
                                 mtcars [, -1])),
                nrow (mtcars))
})

test_that ("LINREG()'s regularisation plots work with a non-default grid", {
  skip_if_not_installed ("glmnet")
  data (mtcars)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  expect_error (LINREG (mtcars [, -1], mtcars [, 1], reg = "lasso",
                        lambda = 10^seq (-3, 3, length.out = 25), graph = TRUE, nrep = 1), NA)
})

# --- Regression: SVR() never read the tuned epsilon back ---------------------------------

test_that ("SVR() reuses the epsilon returned by tune = TRUE", {
  skip_if_not_installed ("e1071")
  data (trees)
  params = SVR (trees [, -3], trees [, 3], tune = TRUE, gamma = 1, cost = 1)
  expect_true ("epsilon" %in% names (params))
  model = SVR (trees [, -3], trees [, 3], methodparameters = params)
  # Reusing pre-tuned parameters must fit a single model, not re-tune from scratch: with
  # epsilon ignored, the length-3 default vector kept triggering a full search.
  expect_equal (model$model$epsilon, params$epsilon)
})

# --- Regression: PCR/PLSR chose their number of components by leave-one-out ---------------

test_that ("LINREG() exposes the validation scheme for PCR and PLSR", {
  skip_if_not_installed ("pls")
  data (mtcars)
  expect_equal (eval (formals (LINREG)$validation) [1], "CV")
  set.seed (1)
  for (r in c ("pcr", "plsr"))
  {
    cv  = LINREG (mtcars [, -1], mtcars [, 1], reg = r, graph = FALSE)
    loo = LINREG (mtcars [, -1], mtcars [, 1], reg = r, validation = "LOO", graph = FALSE)
    expect_true (is.finite (cv$model$ncomp), info = r)
    expect_true (is.finite (loo$model$ncomp), info = r)
    expect_equal (length (predict (cv, mtcars [, -1])), nrow (mtcars), info = r)
  }
})

test_that ("cookplot() and leverageplot() still work without the discarded rstudent() call", {
  data (trees)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  model = LINREG (trees [, -3], trees [, 3], graph = FALSE)
  expect_error (cookplot (model), NA)
  expect_error (leverageplot (model), NA)
  expect_error (resplot (model), NA) # this one genuinely uses rstudent()
})

# =========================================================================================
# Tenth batch: gradient boosting for regression
# =========================================================================================

test_that ("GBREG() fits and predicts a numeric target", {
  skip_if_not_installed ("xgboost")
  data (trees)
  d = splitdata (trees, 3, seed = 0)
  model = GBREG (d$train.x, d$train.y, ntree = 30)
  expect_s3_class (model, "model")
  expect_equal (model$method, "XGBREG")
  pred = predict (model, d$test.x)
  expect_true (is.numeric (pred))
  expect_length (pred, nrow (d$test.x))
  # It is a regression: the predictions live on the scale of the target.
  expect_true (all (pred > 0.5 * min (trees [, 3]) & pred < 2 * max (trees [, 3])))
  expect_gt (evaluation (pred, d$test.y), 0.5)
})

test_that ("GBREG() refuses a categorical target and says which function to use", {
  data (iris)
  expect_error (GBREG (iris [, -5], iris [, 5]), "GRADIENTBOOSTING")
})

test_that ("GBREG() goes through performance() like the other regression methods", {
  skip_if_not_installed ("xgboost")
  data (trees)
  d = splitdata (trees, 3, seed = 0)
  expect_error (performance (GBREG, d$train.x, d$train.y, d$test.x, d$test.y, ntree = 30), NA)
  expect_s3_class (GBREG (d$train.x, d$train.y, tune = TRUE), "params")
})

# =========================================================================================
# Second audit, B2 / B3 / B5
# =========================================================================================

test_that ("the adjusted R2 counts the intercept among the parameters it penalizes", {
  # 1 - (1 - R2) (n - 1) / (n - p - 1), not / (n - p): the intercept is a parameter too. The
  # missing -1 made every value slightly optimistic, and made the penalty vanish exactly where
  # it should bite hardest.
  data (trees)
  model = LINREG (trees [, -3], trees [, 3])
  pred = predict (model, trees [, -3])
  expect_equal (evaluation.adjr2 (pred, trees [, 3], ncol = 2),
                summary (model$model)$adj.r.squared)
  # The three regression criteria and their internal counterparts are the same function, so
  # they cannot drift apart again.
  expect_equal (evaluation.r2 (pred, trees [, 3]), summary (model$model)$r.squared)
  expect_equal (evaluation.msep (pred, trees [, 3]), mean ((pred - trees [, 3])^2))
})

test_that ("the adjusted R2 says what it needs instead of failing on a missing argument", {
  data (trees)
  model = LINREG (trees [, -3], trees [, 3])
  pred = predict (model, trees [, -3])
  expect_error (evaluation (pred, trees [, 3], eval = "adjr2"), "ncol")
  # No residual degree of freedom left: undefined rather than a plausible-looking number.
  expect_warning (value <- evaluation.adjr2 (pred [1:3], trees [1:3, 3], ncol = 2),
                  "undefined")
  expect_true (is.na (value))
  # Through performance(), 'ncol' is supplied automatically.
  d = splitdata (trees, 3, seed = 0)
  expect_true (is.finite (performance (LINREG, d$train.x, d$train.y, d$test.x, d$test.y,
                                       eval = "adjr2")))
})

test_that ("the penalized regressions are fitted under the penalty they selected", {
  skip_if_not_installed ("glmnet")
  data (mtcars)
  x = as.matrix (mtcars [, -1])
  y = mtcars [, 1]
  # lambda was chosen by cross-validation on standardized predictors and the model then
  # fitted on raw ones, so the returned coefficients did not belong to the retained lambda.
  for (reg in c ("ridge", "lasso", "elastic"))
  {
    model = suppressWarnings (LINREG (mtcars [, -1], y, reg = reg, nrep = 2, seed = 0))
    alpha = switch (reg, ridge = 0, lasso = 1, .5)
    reference = glmnet::glmnet (x, y, alpha = alpha, lambda = model$model$lambda,
                                standardize = TRUE)
    expect_equal (as.vector (coef (model$model)), as.vector (coef (reference)), info = reg)
  }
})

test_that ("subset selection defaults to a criterion that can actually reject a variable", {
  skip_if_not_installed ("leaps")
  data (mtcars)
  nvar = function (...)
    length (coef (LINREG (mtcars [, -1], mtcars [, 1], reg = "subset", ...)$model)) - 1
  # The R2 grows with every variable added, so it always keeps them all -- it used to be the
  # default, which made 'subset' select nothing.
  expect_equal (nvar (regeval = "r2"), ncol (mtcars) - 1)
  for (criterion in c ("bic", "adjr2", "cp"))
    expect_lt (nvar (regeval = criterion), ncol (mtcars) - 1)
  expect_equal (nvar (), nvar (regeval = "bic"))
  # pcr and plsr keep their own default, and their own pair of criteria.
  expect_error (LINREG (mtcars [, -1], mtcars [, 1], reg = "pcr"), NA)
  expect_error (LINREG (mtcars [, -1], mtcars [, 1], reg = "plsr", regeval = "msep"), NA)
})

test_that ("subset selection refuses a criterion it cannot compute, by name", {
  skip_if_not_installed ("leaps")
  data (mtcars)
  # 'msep' is offered by the same argument for pcr/plsr; here it used to reach summary()[[NA]]
  # and fail on "subscript out of bounds".
  expect_error (LINREG (mtcars [, -1], mtcars [, 1], reg = "subset", regeval = "msep"),
                "bic, adjr2, cp, r2")
})

# =========================================================================================
# Fourth audit, D1
# =========================================================================================

test_that ("every LINREG algorithm returns the same kind of object", {
  skip_if_not_installed ("leaps")
  skip_if_not_installed ("glmnet")
  data (mtcars)
  # reg = "subset" returned the bare 'lm', so it was the one branch whose result had another
  # class, printed differently, and bypassed predict.model().
  for (reg in c ("linear", "subset", "ridge", "lasso", "elastic", "pcr", "plsr"))
  {
    model = suppressWarnings (LINREG (mtcars [, -1], mtcars [, 1], reg = reg, nrep = 2,
                                      seed = 0))
    expect_s3_class (model, "model", exact = TRUE)
    expect_length (predict (model, mtcars [, -1]), nrow (mtcars))
    expect_output (print (model))
  }
})

test_that ("LINREG names the algorithms and the codings it accepts", {
  data (ToothGrowth)
  # A misspelt 'quali' left the formula unassigned, so paste() picked up stats::formula -- the
  # function -- and the call died on "cannot coerce type 'closure'".
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "intercepte"),
                "'arg' should be one of|should be one of")
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], reg = "ridgge"),
                "'arg' should be one of|should be one of")
  # Partial matching still works, as everywhere else in the package.
  expect_error (LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "inter"), NA)
})

# =========================================================================================
# Sixth audit: MLPREG's output layer
# =========================================================================================

test_that ("MLPREG fits a regression, not a saturating logistic unit", {
  skip_if_not_installed ("nnet")
  data (trees)
  # Without linout = TRUE, nnet() puts a logistic unit on the output. It saturates, and the
  # network predicts very nearly a constant: an R2 of 0.13 on 'trees' where a linear output
  # gives 0.88. The fallback branch of MLPREG asked for it; the two that actually ran did not.
  model = suppressWarnings (MLPREG (trees [, -3], trees [, 3], size = 2, decay = .1, seed = 0))
  pred = predict (model, trees [, -3])
  expect_gt (evaluation.r2 (pred, trees [, 3]), 0.75)
  # A saturated network answers the same thing everywhere; a fitted one spans the target.
  expect_gt (diff (range (pred)), 0.5 * diff (range (trees [, 3])))
  # And the grid search, which goes through tune.nnet(), gets the same treatment.
  tuned = suppressWarnings (MLPREG (trees [, -3], trees [, 3], nfolds = 3, seed = 0))
  expect_gt (evaluation.r2 (predict (tuned, trees [, -3]), trees [, 3]), 0.75)
})

test_that ("the penalty of a regularized regression is chosen once by default", {
  skip_if_not_installed ("glmnet")
  # nrep repetitions of a 10-fold cv.glmnet(), which already averages over its own folds.
  expect_equal (formals (LINREG)$nrep, 1)
  data (mtcars)
  expect_s3_class (suppressWarnings (LINREG (mtcars [, -1], mtcars [, 1], reg = "ridge",
                                             seed = 0)), "model")
})

# --- Ninth audit: the number of components a cross-validation can estimate ------------------

test_that ("PCR and PLS ask for no more components than the validation can fit", {
  skip_if_not_installed ("pls")
  data (cookies.desc.train)
  data (cookies.y.train)
  # 40 observations and 700 variables: asking for nrow - 2 = 38 components made pls warn
  # "`ncomp' reduced to 35 due to cross-validation" on every call -- and p >> n is the very
  # situation these two methods are for.
  for (reg in c ("pcr", "plsr"))
    expect_warning (LINREG (cookies.desc.train, cookies.y.train [, 1], reg = reg), NA,
                    info = reg)
  # Leave-one-out leaves one observation out, so it keeps the old bound.
  expect_warning (LINREG (cookies.desc.train, cookies.y.train [, 1], reg = "pcr",
                          validation = "LOO"), NA)
  # And the model that comes out is a usable one.
  model = LINREG (cookies.desc.train, cookies.y.train [, 1], reg = "pcr")
  data (cookies.desc.test)
  data (cookies.y.test)
  expect_gt (evaluation.r2 (predict (model, cookies.desc.test), cookies.y.test [, 1]), 0.8)
  # A small dataset is unaffected.
  data (trees)
  expect_warning (LINREG (trees [, -3], trees [, 3], reg = "plsr"), NA)
})
