#' Plot the Cook's distance of a linear regression model
#'
#' Plot the Cook's distance of a linear regression model.
#'
#' @name cookplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for the x-axis.
#' @param labels The labels of the instances.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' cookplot (model)
cookplot <-
  function (model, index = NULL, labels = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    if (is.null (labels))
      labels = 1:length (mod$residuals)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else if (length (index) == 1)
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    else
      xlab = ""
    y = stats::cooks.distance (mod)
    n = length (mod$residuals)
    threshold = 4 / n
    ylim = c (0, max (c (y, threshold)) * 1.11)
    graphics::plot (index, y, type = "h", ylim = ylim, ylab = "Cook's distance", xlab = xlab, cex = 2, lwd = 2, col = ifelse (y <= threshold, 1, 2), cex.axis = 1.5, cex.lab = 1.5)
    graphics::abline (h = threshold, lty = 2, lwd = 2)
    select = which (y > threshold)
    if (length (select) > 0)
      graphics::text (index [select], y [select], labels [select], pos = 3, cex = 1)
  }

#' @keywords internal
# The adjusted R2 penalises the R2 by the number of *parameters* the model spends, which is
# one per predictor plus the intercept: the residual degrees of freedom are n - p - 1, not
# n - p. The missing intercept made every value slightly too optimistic (0.9462 instead of
# summary (lm)'s 0.9442 on 'trees') and, worse, made the penalty vanish exactly when it should
# bite hardest -- with as many predictors as observations it read 1 - 0 = 1 rather than
# undefined.
eval.adjr2 <-
  function (predictions, gt, nrow = length (predictions), ncol, ...)
  {
    if (missing (ncol) || is.null (ncol))
      stop ("evaluation.adjr2: 'ncol' (the number of predictors the model uses) is needed to ",
            "penalize the R2 and has no default. Pass it explicitly, e.g. ",
            "evaluation.adjr2 (pred, gt, ncol = ncol (test.x)); through evaluation() or ",
            "performance() it is supplied for you.")
    df = nrow - ncol - 1
    if (df <= 0)
    {
      warning ("evaluation.adjr2: the model uses ", ncol, " predictor(s) for ", nrow,
               " observation(s), leaving ", df, " residual degree(s) of freedom: the adjusted ",
               "R2 is undefined.")
      return (NA_real_)
    }
    return (1 - (1 - eval.r2 (predictions, gt)) * (nrow - 1) / df)
  }

#' @keywords internal
eval.msep <-
  function (predictions, gt, ...) mean ((predictions - gt)^2)

#' @keywords internal
eval.r2 <-
  function (predictions, gt, ...) 1 - ((sum ((predictions - gt)^2) / sum ((gt - mean (gt))^2)))

#' Adjusted R2 evaluation of regression predictions
#'
#' Evaluation predictions of a regression model according to the adjusted R2, i.e. the R2
#' penalized by the number of variables used by the model.
#' @name evaluation.adjr2
#' @param predictions The predictions of a regression model (\code{vector}).
#' @param gt The ground truth (\code{vector}).
#' @param nrow Number of observations (defaults to the number of predictions).
#' @param ncol Number of predictors used by the model. This one has no default: the adjustment
#' cannot be computed without it. The residual degrees of freedom are \code{nrow - ncol - 1},
#' the intercept counting as one parameter; the value is \code{NA} when they run out.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.r2}}, \code{\link{evaluation.msep}}, \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (trees)
#' d = splitdata (trees, 3)
#' model.linreg = LINREG (d$train.x, d$train.y)
#' pred.linreg = predict (model.linreg, d$test.x)
#' evaluation.adjr2 (pred.linreg, d$test.y, ncol = ncol (d$test.x))
evaluation.adjr2 <-
  function (predictions, gt, nrow = length (predictions), ncol, ...)
    eval.adjr2 (predictions, gt, nrow, ncol)

#' MSEP evaluation of regression predictions
#'
#' Evaluation predictions of a regression model according to MSEP
#' @name evaluation.msep
#' @param predictions The predictions of a regression model (\code{vector}).
#' @param gt The ground truth (\code{vector}).
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.r2}}, \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (trees)
#' d = splitdata (trees, 3)
#' model.lin = LINREG (d$train.x, d$train.y)
#' pred.lin = predict (model.lin, d$test.x)
#' evaluation.msep (pred.lin, d$test.y)
evaluation.msep <-
  function (predictions, gt, ...) eval.msep (predictions, gt)

#' R2 evaluation of regression predictions
#'
#' Evaluation predictions of a regression model according to R2
#' @name evaluation.r2
#' @param predictions The predictions of a regression model (\code{vector}).
#' @param gt The ground truth (\code{vector}).
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.msep}}, \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (trees)
#' d = splitdata (trees, 3)
#' model.linreg = LINREG (d$train.x, d$train.y)
#' pred.linreg = predict (model.linreg, d$test.x)
#' evaluation.r2 (pred.linreg, d$test.y)
evaluation.r2 <-
  function (predictions, gt, ...) eval.r2 (predictions, gt)

#' Kernel Regression
#'
#' This function builds a kernel regression model.
#' @name KERREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param bandwidth The bandwidth parameter.
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{KERREG} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{KERREG} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[ibr]{npregress}}
#' @examples
#' require (datasets)
#' data (trees)
#' KERREG (trees [, -3], trees [, 3])
KERREG <-
  function (x, y, bandwidth = 1,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      res = list (model = ibr::npregress (y = y, x, bandwidth = bandwidth), method = "KERREG")
      class (res) = "model"
    }
    return (res)
  }


#' Plot the leverage points of a linear regression model
#'
#' Plot the leverage points of a linear regression model.
#'
#' @name leverageplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for the x-axis.
#' @param labels The labels of the instances.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' leverageplot (model)
leverageplot <-
  function (model, index = NULL, labels = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    if (is.null (labels))
      labels = 1:length (mod$residuals)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else if (length (index) == 1)
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    else
      xlab = ""
    y = stats::hatvalues (mod)
    p = length (mod$coefficients)
    n = length (mod$residuals)
    thresholds = c (2 * p / n, 3 * p / n)
    ylim = c (0, max (c (y, thresholds)) * 1.11)
    graphics::plot (index, y, type = "h", ylim = ylim, ylab = expression('h'['oo']), xlab = xlab, cex = 2, lwd = 2, col = ifelse (y <= thresholds [1], 1, 2), cex.axis = 1.5, cex.lab = 1.5)
    graphics::abline (h = thresholds, lty = 2:3, lwd = 2)
    select = which (y > thresholds [1])
    if (length (select) > 0)
      graphics::text (index [select], y [select], labels [select], pos = 3, cex = 1)
  }

#' Linear Regression
#'
#' This function builds a linear regression model.
#' Standard least square method, variable selection, factorial methods are available.
#' @name LINREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param quali Indicates how to use the qualitative variables.
#' @param reg The algorithm.
#' @param regeval The criterion used to choose between models. For \code{reg = "subset"}:
#' \code{"bic"} (the default), \code{"adjr2"} or \code{"cp"}, which all penalize the number of
#' variables, and \code{"r2"}, which does not -- the R2 can only grow when a variable is added,
#' so it always retains every variable. For \code{reg = "pcr"} and \code{reg = "plsr"}, where
#' the choice is a number of components: \code{"r2"} (the default) or \code{"msep"}. Ignored by
#' the other algorithms.
#' @param scale If true, PCR and PLS use scaled dataset.
#' @param validation How the number of components of a PCR or PLS regression is chosen:
#' \code{"CV"} (10 random segments, the default) or \code{"LOO"} (leave-one-out, one
#' regression per observation -- only practical on a small dataset). Ignored by the other
#' algorithms.
#' @param lambda The lambda parameter of Ridge, Lasso and Elastic net regression.
#' @param alpha The elasticnet mixing parameter.
#' @param nrep How many times the cross-validation choosing \code{lambda} is repeated, its
#' errors being averaged. One is enough in practice -- \code{\link[glmnet]{cv.glmnet}} already
#' averages over its own folds -- and each extra repetition costs another ten fits. Raise it to
#' steady the choice of \code{lambda} on a small or noisy dataset.
#' @param graph A logical indicating whether or not graphics should be plotted (ridge, LASSO and elastic net).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{LINREG} does not
#' support reusing pre-tuned parameters.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[stats]{lm}}, \code{\link[leaps]{regsubsets}}, \code{\link[pls]{mvr}}, \code{\link[glmnet]{glmnet}}
#' @examples
#' \dontrun{
#' require (datasets)
#' # With one independent variable
#' data (cars)
#' LINREG (cars [, -2], cars [, 2])
#' # With two independent variables
#' data (trees)
#' LINREG (trees [, -3], trees [, 3])
#' # With non numeric variables
#' data (ToothGrowth)
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "intercept") # Different intercept
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "slope") # Different slope
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], quali = "both") # Complete model
#' # With multiple numeric variables
#' data (mtcars)
#' LINREG (mtcars [, -1], mtcars [, 1])
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "subset", regeval = "adjr2")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "ridge")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "lasso")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "elastic")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "pcr")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "plsr")
#' }
LINREG <-
  function (x, y, quali = c ("none", "intercept", "slope", "both"),
            reg = c ("linear", "subset", "ridge", "lasso", "elastic", "pcr", "plsr"),
            # The two algorithms that use a criterion do not accept the same ones, and do not
            # want the same default. Subset selection compares models of *different sizes*, so
            # its criterion has to penalize size: the R2 cannot, since it only grows when a
            # variable is added, and using it selected all 11 predictors of 'mtcars' every time
            # -- a variable selection that never selects anything.
            regeval = if (reg [1] == "subset") c ("bic", "adjr2", "cp", "r2")
                      else c ("r2", "msep"),
            scale = TRUE,
            validation = c ("CV", "LOO"),
            lambda = 10^seq (-5, 5, length.out = 101),
            alpha = .5,
            nrep = 1,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    model = NULL
    if (tune)
      model = emptyparams ()
    else
    {
      if (is.vector (x))
        x = data.frame (X = x)
      reg = match.arg (reg [1], c ("linear", "subset", "ridge", "lasso", "elastic", "pcr", "plsr"))
      quali = match.arg (quali [1], c ("none", "intercept", "slope", "both"))
      model = NULL
      if (reg [1] == "linear")
      {
        isquali = sapply (as.data.frame (x), is.factor)
        if (!any (isquali))
          quali = "none"
        nquali = colnames (x) [isquali]
        if (sum (isquali) > 0)
          nquali = paste ("`", nquali, "`", sep = "")
        nquanti = colnames (x) [!isquali]
        if (sum (!isquali) > 0)
          nquanti = paste ("`", nquanti, "`", sep = "")
        lquali = paste (nquali, collapse = "+")
        lquanti = paste (nquanti, collapse = "+")
        if (quali [1] == "none")
          formula = lquanti
        if (quali [1] == "intercept")
          formula = paste ("-1", lquali, lquanti, sep = "+")
        if (quali [1] == "slope")
          formula = paste (apply (expand.grid (nquali, nquanti), 1, paste, collapse = ":"), collapse = "+")
        if (quali [1] == "both")
          formula = paste ("-1", lquali, paste (apply (expand.grid (nquali, nquanti), 1, paste, collapse = ":"), collapse = "+"), sep = "+")
        f = stats::as.formula (paste ("y~", formula, sep = ""))
        model = stats::lm (formula = f, x)
        model = list (model = model, method = "lm")
        class (model) = "model"
      }
      else if (reg [1] == "subset")
      {
        # Wrapped like every other branch, so that every value of 'reg' returns the same kind
        # of object.
        model = list (model = lmselect (x, y, regeval, graph, ...), method = "lm")
        class (model) = "model"
      }
      else if ((reg [1] == "ridge") | (reg [1] == "lasso") | (reg [1] == "elastic"))
      {
        palpha = alpha
        if (reg [1] == "ridge")
          palpha = 0
        else if (reg [1] == "lasso")
          palpha = 1
        # glmnet always works on a decreasing lambda sequence, and that is the order of the
        # cross-validation errors it returns.
        lambda = sort (lambda, decreasing = TRUE)
        cv = NULL
        for (i in 1:nrep)
          cv = rbind (cv, glmnet::cv.glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda, standardize = TRUE, nfolds = 10)$cvm)
        cv = apply (cv, 2, mean)
        lambda.min = lambda [which.min (cv)]
        # standardize = TRUE everywhere, as in the cross-validation just above: the penalty
        # lambda.min stands for has to be the penalty the returned model is fitted under.
        # glmnet reports the coefficients on the original scale either way; 'standardize' only
        # says whether each predictor is penalized in its own units or in standard deviations,
        # and the latter is the only reading under which one lambda means the same thing for
        # all of them.
        if (graph)
        {
          plotmsep (lambda, cv)
          plotcoeff (lambda, glmnet::glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda, standardize = TRUE), cv)
        }
        model = list (model = glmnet::glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda.min, standardize = TRUE), method = "regularization")
        class (model) = "model"
      }
      else if ((reg [1] == "pcr") | (reg [1] == "plsr"))
      {
        model = list (model = lmfact (x = x, y = y, reg = reg [1], regeval = regeval, scale = scale, graph = graph, validation = validation [1], ...), method = "MRV")
        class (model) = "model"
      }
      else
        stop ("Invalid regression method")
    }
    return (model)
  }

#' @keywords internal
lmfact <-
  function (x, y, reg, regeval = c ("r2", "msep"), scale = TRUE, ...)
  {
    nbc = nbcomp (x = x, y = y, reg = reg, regeval = regeval, ...)
    # switch() rather than get (reg): the same named-lookup treatment the rest of the package
    # received, and nbcomp() just above already had it.
    regfun = switch (reg [1], pcr = pls::pcr, plsr = pls::plsr,
                     stop ("lmfact: unsupported 'reg' method: ", reg [1]))
    model = regfun (y~., data = x, ncomp = nbc, scale = scale)
    return (model)
  }

#' @keywords internal
# Each criterion goes under three names: the one the package uses, the one plot.regsubsets()
# wants for its 'scale' argument, and the one summary.regsubsets() gives its component.
# 'sign' is +1 for a criterion to maximise, -1 for one to minimise.
lmselect.criteria <-
  function ()
    list (bic   = list (scale = "bic",   field = "bic",   sign = -1),
          adjr2 = list (scale = "adjr2", field = "adjr2", sign =  1),
          cp    = list (scale = "Cp",    field = "cp",    sign = -1),
          r2    = list (scale = "r2",    field = "rsq",   sign =  1))

#' @keywords internal
lmselect <-
  function (x, y, regeval = c ("bic", "adjr2", "cp", "r2"), graph = TRUE, ...)
  {
    criteria = lmselect.criteria ()
    name = tolower (regeval [1])
    if (name == "rsq")
      name = "r2"
    crit = criteria [[name]]
    # 'msep' is offered by the same argument for pcr and plsr, but regsubsets() cannot
    # compute it.
    if (is.null (crit))
      stop ("LINREG: reg = \"subset\" cannot use regeval = \"", regeval [1],
            "\". Available criteria: ", paste (names (criteria), collapse = ", "), ".")
    rss = leaps::regsubsets (y~., x, nvmax = ncol (x), method = "exhaustive")
    if (graph)
      graphics::plot (rss, scale = crit$scale)
    s = summary (rss)
    best = which.max (crit$sign * s [[crit$field]])
    var = colnames (s$which) [s$which [best, ]][-1]
    f = stats::reformulate (var, response = "y")
    model = stats::lm (f, data = x)
    return (model)
  }

#' Multi-Layer Perceptron Regression
#'
#' This function builds a regression model using MLP.
#' @name MLPREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param size The size of the hidden layer (if a vector, cross-over validation is used to chose the best size).
#' @param decay The decay (between 0 and 1) of the backpropagation algorithm (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{size} and
#' \code{decay}. Named (and behaves identically to) \code{methodparameters} rather than
#' \code{params}, for consistency with \code{\link{MLP}} and with the calling convention used by
#' \code{\link{performance}}/the internal \code{protocol.*} functions, which always pass a
#' \code{methodparameters} argument.
#' @inheritParams tune.doc
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{MLPREG} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[nnet]{nnet}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (trees)
#' MLPREG (trees [, -3], trees [, 3])
#' }
MLPREG <-
  function (x,
            y,
            # if/else rather than 2:(ifelse (...)), to read the same way as MLP()'s 'hidden'
            # default and to keep ncol (x) from being evaluated when x is a vector. Same values
            # as before.
            size = if (is.vector (x)) 2 else 2:ncol (x),
            decay = 10^(-3:-1),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    # No early 'return (emptyparams ())' here: tune = TRUE must return the size and decay this
    # method actually retained, exactly as MLP() does, so that performance() can obtain them
    # once and hand them back on every split instead of re-tuning the network each time. The
    # params object is built at the end of the body, once the search has run.
    setseed (seed)
    model = NULL
    d = cbind.data.frame (Y = y, x)
    minimum = apply (d, 2, min)
    range = apply (d, 2, max) - apply (d, 2, min)
    range [range == 0] = 1 # Avoid division by zero for constant variables (NaN/Inf)
    d.norm = sweep (sweep (d, 2, minimum, FUN = "-"), 2, range, FUN = "/")
    if (is.null (colnames (x)))
      colnames (d.norm) = c ("Y", paste ("X", 1:(ncol (d.norm) - 1), sep = ""))
    else
      colnames (d.norm) = c ("Y", colnames (x))
    if (!is.null (methodparameters))
    {
      # Each value is taken only if the 'params' object actually carries it: an *empty* one
      # -- what a method with nothing to tune returns, and what performance() then hands back
      # to it -- must leave the defaults alone rather than overwrite them with NULL.
      if (!is.null (methodparameters$hidden))
        size = methodparameters$hidden
      if (!is.null (methodparameters$decay))
        decay = methodparameters$decay
    }
    if (length (size) > 1 | length (decay) > 1)
    {
      # linout = TRUE: this is a regression, so the output unit is linear. Without it nnet()
      # puts a logistic one there, which saturates -- the network then predicts very nearly a
      # constant. Fitting one hidden unit too few on 'trees' that way gives an R2 of 0.13 where
      # a linear output gives 0.88. The fallback below already asked for it; the two branches
      # that actually run did not.
      tuned = e1071::tune.nnet (Y~., data = d.norm, size = size, decay = decay,
                                tunecontrol = tune.scheme (nfolds), linout = TRUE, ...)
      model = tuned$best.model
      # See SVM(): tune.*() leaves 'best.model' empty when it could not refit one, and the
      # NULL then travelled to predict().
      if (is.null (model))
        model = nnet::nnet (Y~., data = d.norm, size = tuned$best.parameters$size,
                            decay = tuned$best.parameters$decay, linout = TRUE, trace = FALSE)
    }
    else
      model = nnet::nnet (Y~., data = d.norm, size = size, decay = decay, linout = TRUE,
                          trace = FALSE, ...)
    res = NULL
    if (tune)
    {
      res = list (decay = model$decay, hidden = model$n [2])
      class (res) = "params"
    }
    else
    {
      res = list (model = model, minimum = minimum, range = range)
      res = list (model = res, method = "MLPREG")
      class (res) = "model"
    }
    return (res)
  }

#' @keywords internal
nbcomp <-
  function (x, y, reg, regeval = c ("r2", "msep"), graph = TRUE, validation = "CV", ...)
  {
    eval = toupper (regeval [1])
    optim = 0
    if (eval [1] == "R2")
      optim = 1
    else if (eval [1] == "MSEP")
      optim = -1
    regfun = switch (reg [1], pcr = pls::pcr, plsr = pls::plsr, stop ("nbcomp: unsupported 'reg' method: ", reg [1]))
    evalfun = switch (eval [1], R2 = pls::R2, MSEP = pls::MSEP, stop ("nbcomp: unsupported evaluation criterion: ", eval [1]))
    # "CV" (10 random segments, the pls default) rather than "LOO", which fits one regression
    # per observation: both the usual choice and orders of magnitude cheaper on anything but a
    # small dataset. LINREG's 'validation' argument asks for "LOO" when that is what is wanted.
    # A cross-validation fits the model on *segments* of the data, so the number of
    # components it can estimate is bounded by the size of a training segment, not by the whole
    # sample. Asking for nrow (x) - 2 made pls warn ("`ncomp' reduced to 35 due to
    # cross-validation") on every call with more variables than observations -- which is the
    # very situation PCR and PLS are for.
    n = nrow (x)
    segments = if (validation [1] == "LOO") n else min (10, n)
    largest = n - ceiling (n / segments) - 1
    ncomp = max (1, min (ncol (x), n - 2, largest))
    fit = regfun (y~., data = x, ncomp = ncomp, validation = validation)
    res = evalfun (fit, estimate = c ("train", "CV"))
    if (optim < 0)
      ncomp = which.min (res$val ["CV",,]) - 1
    else
      ncomp = which.max (res$val ["CV",,]) - 1
    ncomp = max (1, as.numeric (ncomp))
    if (graph)
    {
      xlab = paste ("Number of components (", substr (toupper (reg), 1, 3), ")", sep = "");
      graphics::plot (res, main = "", ylab = eval, xlab = xlab)
      pos = "topright"
      if (optim > 0)
        pos = "bottomright"
      graphics::legend (pos, c ("Learning set",
                                if (validation [1] == "LOO") "LOOCV" else "Cross-validation"),
                        lty = 1:2, col = 1:2)
      graphics::abline (v = ncomp, lty = 2)
      bottom = min (res$val ["CV",,])
      top = max (res$val ["CV",,])
      delta = top - bottom
      graphics::text (ncomp, bottom + delta * .5, paste ("Nb. comp. :", ncomp), pos = 4)
    }
    return (ncomp)
  }

#' Plot actual vs. predictions
#'
#' Plot actual vs. predictions of a regression model.
#' @name plotavsp
#' @param predictions The predictions of a classification model (\code{vector}).
#' @param gt The ground truth of the dataset (\code{vector}).
#' @export
#' @seealso \code{\link{confusion}}, \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation.msep}}, \code{\link{evaluation.r2}}, \code{\link{performance}}
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' pred = predict (model, trees [, -3])
#' plotavsp (pred, trees [, 3])
plotavsp <-
  function (predictions, gt)
  {
    palette = grDevices::colorRampPalette (c ("forestgreen", "red")) (101)
    sqres = (gt - predictions)^2 / mean ((gt - mean (gt))^2)
    col = palette [1 + round (100 * sapply (sqres, FUN = function (x) min (x, 1)), 0)]
    plot (gt, predictions, asp = 1, xlab = "Actual", ylab = "Predicted", col = col)
    graphics::abline (0, 1, col = "forestgreen")
  }

#' @keywords internal
plotcoeff <-
  function (lambda, model, cv = NULL)
  {
    # 'lambda' is already sorted decreasing by its caller, matching the column order of
    # model$beta and of the cross-validation errors.
    x = log (lambda)
    y = t (model$beta)
    graphics::matplot (x, y, type = "l", lwd = 2, xlab = expression (paste ("log(", lambda, ")")), ylab = "Coefficients",
                       cex.lab = 1.5, cex.axis = 1.5, col = 1:9, lty = 1:9)
    bestx = x [which.min (cv)]
    graphics::abline (v = bestx, lwd = 2, lty = 2, col = "darkgrey")
    graphics::legend ("topright", colnames (y), col = 1:9, lty = 1:9, lwd = 2)
  }

#' @keywords internal
plotmsep <-
  function (lambda, cv)
  {
    x = log (lambda, base = 10)
    y = cv
    graphics::plot (x, y, t = "l", lwd = 2, xlab = expression(paste("log(", lambda, ")")), ylab = "MSEP",
                    cex.lab = 1.5, cex.axis = 1.5, col = "red")
    bestx = x [which.min (y)]
    graphics::abline (v = bestx, lwd = 2, lty = 2, col = "darkgrey")
    graphics::text (bestx, max (y) * .95, bquote(paste("log(",lambda, ")=", .(bestx))), pos = 2, cex = 1.5)
    graphics::text (bestx, max (y) * .85, bquote(paste(lambda, "=", .(lambda [which.min (y)]))), pos = 2, cex = 1.5)
  }

#' Polynomial Regression
#'
#' This function builds a polynomial regression model.
#' @name POLYREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param degree The polynom degree.
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{POLYREG} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{POLYREG} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[mda]{polyreg}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (trees)
#' POLYREG (trees [, -3], trees [, 3])
#' }
POLYREG <-
  function (x, y, degree = 2,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      res = list (model = mda::polyreg (x, y, degree = degree), method = "POLYREG")
      class (res) = "model"
    }
    return (res)
  }

#' Plot function for a regression model
#'
#' Plot a regression model on a 2-D plot. The predictor \code{x} should be one-dimensional.
#'
#' @name regplot
#' @param model The model to be plotted.
#' @param x The predictor \code{vector}.
#' @param y The response \code{vector}.
#' @param margin A margin parameter.
#' @param ... Other graphical parameters
#' @export
#' @examples
#' require (datasets)
#' data (cars)
#' model = POLYREG (cars [, -2], cars [, 2])
#' regplot (model, cars [, -2], cars [, 2])
regplot <-
  function (model, x, y, margin = .1, ...)
  {
    mod = model
    if (methods::is (mod, "model"))
      mod = mod$model
    deltax = (max (x) - min (x)) * margin
    xlim = c (min (x) - deltax, max (x) + deltax)
    deltay = (max (y) - min (y)) * margin
    ylim = c (min (y) - deltay, max (y) + deltay)
    xl = data.frame (X = seq (xlim [1], xlim [2], length = 1000))
    if (!is.null (mod$terms))
      colnames (xl) = attr (attr (mod$terms, "factor"), "dimnames") [[1]] [2]
    graphics::plot (x, y, xaxs = "i", xlim = xlim, ylim = ylim, ...)
    graphics::lines (cbind (xl, stats::predict (model, xl)), col = 2)
  }

#' Plot the studentized residuals of a linear regression model
#'
#' Plot the studentized residuals of a linear regression model.
#'
#' @name resplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for the x-axis.
#' @param labels The labels of the instances.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' resplot (model) # Ordered by index
#' resplot (model, index = 0) # Ordered by variable "Volume" (dependent variable)
#' resplot (model, index = 1) # Ordered by variable "Girth" (independent variable)
#' resplot (model, index = 2) # Ordered by variable "Height" (independent variable)
resplot <-
  function (model, index = NULL, labels = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    y = stats::rstudent (mod)
    if (is.null (labels))
      labels = 1:length (mod$residuals)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else if (length (index) == 1)
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    else
      xlab = ""
    ylim = c (min (y, -2) * 1.11, max (y, 2) * 1.11)
    graphics::plot (index, y, ylim = ylim, ylab = "Residuals", xlab = xlab, cex = 2, lwd = 2, col = ifelse (abs (y) <= 2, "darkgray", 2), cex.axis = 1.5, cex.lab = 1.5)
    mod = stats::loess (y ~ index)
    x = seq (min (index), max (index), length.out = 1001)
    graphics::lines (x, stats::predict (mod, x), lwd = 2, lty = 4, col = 4)
    graphics::abline (h = c (-2, 0, 2), lty = c (2, 1, 2), lwd = 2)
    select = which (abs (y) > 2)
    if (length (select) > 0)
      graphics::text (index [select], y [select], labels = labels [select], pos = 3, cex = 1)
  }

#' Regression using Support Vector Machine
#'
#' This function builds a regression model using Support Vector Machine.
#' @name SVR
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param gamma The gamma parameter (if a vector, cross-over validation is used to chose the best size).
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param kernel The kernel type.
#' @param epsilon The epsilon parameter (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{epsilon},
#' \code{gamma} and \code{cost}. Named (and behaves identically to) \code{methodparameters} rather
#' than \code{params}, for consistency with \code{\link{SVM}} and with the calling convention used
#' by \code{\link{performance}}/the internal \code{protocol.*} functions, which always pass a
#' \code{methodparameters} argument.
#' @inheritParams tune.doc
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{SVR} does not produce a plot.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVRl}}, \code{\link{SVRr}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (trees)
#' SVR (trees [, -3], trees [, 3], kernel = "linear", cost = 1)
#' SVR (trees [, -3], trees [, 3], kernel = "radial", gamma = 1, cost = 1)
#' }
SVR <-
  function (x,
            y,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            kernel = c ("radial", "linear"),
            epsilon = c (.1, .5, 1),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    check.numeric.predictors (x, "SVR")
    model = NULL
    if (!is.null (methodparameters))
    {
      # Each value is taken only if the 'params' object actually carries it: an *empty* one
      # -- what a method with nothing to tune returns, and what performance() then hands back
      # to it -- must leave the defaults alone rather than overwrite them with NULL.
      if (!is.null (methodparameters$gamma))
        gamma = methodparameters$gamma
      if (!is.null (methodparameters$cost))
        cost = methodparameters$cost
      # epsilon is read back like gamma and cost -- tune = TRUE returns all three. Without it
      # the default length-3 grid keeps the test below TRUE and a full re-tuning runs anyway.
      if (!is.null (methodparameters$epsilon))
        epsilon = methodparameters$epsilon
    }
    if (kernel [1] == "linear")
      gamma = 0
    if (length (gamma) > 1 | length (cost) > 1 | length (epsilon) > 1)
    {
      tuned = e1071::tune.svm (x, y, epsilon = epsilon,
                               gamma = gamma, cost = cost, kernel = kernel [1],
                               tunecontrol = tune.scheme (nfolds), ...)
      model = tuned$best.model
      # See SVM().
      if (is.null (model))
        model = e1071::svm (x, y, epsilon = tuned$best.parameters$epsilon,
                            gamma = tuned$best.parameters$gamma,
                            cost = tuned$best.parameters$cost, kernel = kernel [1], ...)
    }
    else
      model = e1071::svm (x, y, epsilon = epsilon, gamma = gamma, cost = cost, kernel = kernel [1], ...)
    res = NULL
    if (tune)
    {
      res = list (epsilon = model$epsilon, gamma = model$gamma, cost = model$cost)
      class (res) = "params"
    }
    else
    {
      res = list (model = model, method = "SVR")
      class (res) = "model"
    }
    return (res)
  }

#' Regression using Support Vector Machine with a linear kernel
#'
#' This function builds a regression model using Support Vector Machine with a linear kernel.
#' @name SVRl
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param epsilon The epsilon parameter (if a vector, cross-over validation is used to chose the best size).
#' @inheritParams tune.doc
#' @param methodparameters Object containing the parameters. If given, it replaces \code{epsilon}, \code{gamma} and \code{cost}. Named to match \code{\link{SVR}} (see there).
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVR}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (trees)
#' SVRl (trees [, -3], trees [, 3], cost = 1)
#' }
SVRl <-
  function (x,
            y,
            cost = 2^(-3:3),
            epsilon = c (.1, .5, 1),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    return (SVR (
      x = x,
      y = y,
      cost = cost,
      epsilon = epsilon,
      kernel = "linear",
      nfolds = nfolds,
      tune = tune,
      methodparameters = methodparameters,
      graph = graph,
      seed = seed,
      ...
    ))
  }

#' Regression using Support Vector Machine with a radial kernel
#'
#' This function builds a regression model using Support Vector Machine with a radial kernel.
#' @name SVRr
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param gamma The gamma parameter (if a vector, cross-over validation is used to chose the best size).
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param epsilon The epsilon parameter (if a vector, cross-over validation is used to chose the best size).
#' @inheritParams tune.doc
#' @param methodparameters Object containing the parameters. If given, it replaces \code{epsilon}, \code{gamma} and \code{cost}. Named to match \code{\link{SVR}} (see there).
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVR}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (trees)
#' SVRr (trees [, -3], trees [, 3], gamma = 1, cost = 1)
#' }
SVRr <-
  function (x,
            y,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            epsilon = c (.1, .5, 1),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    return (SVR (
      x = x,
      y = y,
      gamma = gamma,
      cost = cost,
      epsilon = epsilon,
      kernel = "radial",
      nfolds = nfolds,
      tune = tune,
      methodparameters = methodparameters,
      graph = graph,
      seed = seed,
      ...
    ))
  }
