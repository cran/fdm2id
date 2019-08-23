#' Plot the Cook's distance of a linear regression model
#'
#' Plot the Cook's distance of a linear regression model.
#'
#' @name cookplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for for the x-axis.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' cookplot (model)
cookplot <-
  function (model, index = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    y = stats::rstudent (mod)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    y = stats::cooks.distance (mod)
    n = length (mod$residuals)
    threshold = 4 / n
    ylim = c (0, max (c (y, threshold)) * 1.11)
    graphics::plot (index, y, type = "h", ylim = ylim, ylab = "Cook's distance", xlab = xlab, cex = 2, lwd = 2, col = ifelse (y <= threshold, 1, 2), cex.axis = 1.5, cex.lab = 1.5)
    graphics::abline (h = threshold, lty = 2, lwd = 2)
    select = which (y > threshold)
    if (length (select) > 0)
      graphics::text (index [select], y [select], select, pos = 3, cex = 1)
  }

#' @keywords internal
eval.msep <-
  function (predictions, targets, ...) mean ((predictions - targets)^2)

#' @keywords internal
eval.r2 <-
  function (predictions, targets, ...) 1 - ((sum ((predictions - targets)^2) / sum ((targets - mean (targets))^2)))

#' MSEP evaluation of regression predictions
#'
#' Evaluation predictions of a regression model according to MSEP
#' @name evaluation.msep
#' @param predictions The predictions of a regression model (\code{vector}).
#' @param targets Actual targets of the dataset (\code{vector}).
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
  function (predictions, targets) mean ((predictions - targets)^2)

#' R2 evaluation of regression predictions
#'
#' Evaluation predictions of a regression model according to R2
#' @name evaluation.r2
#' @param predictions The predictions of a regression model (\code{vector}).
#' @param targets Actual targets of the dataset (\code{vector}).
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
  function (predictions, targets) 1 - ((sum ((predictions - targets)^2) / sum ((targets - mean (targets))^2)))

#' Kernel Regression
#'
#' This function builds a kernel regression model.
#' @name KERREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param bandwidth The bandwidth parameter.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[ibr]{npregress}}
#' @examples
#' require (datasets)
#' data (trees)
#' KERREG (trees [, -3], trees [, 3])
KERREG <-
  function (x, y, bandwidth = 1, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
      res = ibr::npregress (y = y, x, bandwidth = bandwidth)
    return (res)
  }


#' Plot the leverage points of a linear regression model
#'
#' Plot the leverage points of a linear regression model.
#'
#' @name leverageplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for for the x-axis.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' leverageplot (model)
leverageplot <-
  function (model, index = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    y = stats::rstudent (mod)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    y = stats::hatvalues (mod)
    p = length (mod$coefficients)
    n = length (mod$residuals)
    thresholds = c (2 * p / n, 3 * p / n)
    ylim = c (0, max (c (y, thresholds)) * 1.11)
    graphics::plot (index, y, type = "h", ylim = ylim, ylab = expression('h'['oo']), xlab = xlab, cex = 2, lwd = 2, col = ifelse (y <= thresholds [1], 1, 2), cex.axis = 1.5, cex.lab = 1.5)
    graphics::abline (h = thresholds, lty = 2:3, lwd = 2)
    select = which (y > thresholds [1])
    if (length (select) > 0)
      graphics::text (index [select], y [select], select, pos = 3, cex = 1)
  }

#' Linear Regression
#'
#' This function builds a linear regression model.
#' Standard least square method, variable selection, factorial methods are available.
#' @name LINREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param formula A symbolic description of the model to be fitted (as a character string).
#' @param reg The algorithm.
#' @param regeval The evaluation criterion for subset selection.
#' @param scale If true, PCR and PLS use scaled dataset.
#' @param lambda The lambda parameter of Ridge, Lasso and Elastic net regression.
#' @param alpha The elasticnet mixing parameter.
#' @param graph A logical indicating whether or not graphics should be plotted (ridge, LASSO and elastic net).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[stats]{lm}}, \code{\link[leaps]{regsubsets}}, \code{\link[pls]{mvr}}, \code{\link[glmnet]{glmnet}}
#' @examples
#' require (datasets)
#' # With one independant variable
#' data (cars)
#' LINREG (cars [, -2], cars [, 2])
#' # With two independant variables
#' data (trees)
#' LINREG (trees [, -3], trees [, 3])
#' # With non numeric variables
#' data (ToothGrowth)
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], formula = "-1+supp+dose") # Different intersept
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], formula = "dose:supp") # Different slope
#' LINREG (ToothGrowth [, -1], ToothGrowth [, 1], formula = "-1+supp+dose:supp") # Complete model
#' # With multiple numeric variables
#' data (mtcars)
#' LINREG (mtcars [, -1], mtcars [, 1])
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "subset", regeval = "adjr2")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "ridge")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "lasso")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "elastic")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "pcr")
#' LINREG (mtcars [, -1], mtcars [, 1], reg = "plsr")
LINREG <-
  function (x, y, formula = ".",
            reg = c ("linear", "subset", "ridge", "lasso", "elastic", "pcr", "plsr"),
            regeval = c ("r2", "bic", "adjr2", "cp", "msep"),
            scale = TRUE,
            lambda = 10^seq (-5, 5, length.out = 101),
            alpha = .5,
            graph = TRUE,
            tune = FALSE, ...)
  {
    model = NULL
    if (tune)
      model = emptyparams ()
    else
    {
      if (is.vector (x))
        x = data.frame (X = x)
      model = NULL
      if (reg [1] == "linear")
      {
        f = stats::as.formula (paste ("y~", formula, sep = ""))
        model = stats::lm (formula = f, x)
        model = list (model = model, method = "lm")
        class (model) = "model"
      }
      else if (reg [1] == "subset")
        model = lmselect (x, y, regeval, graph, ...)
      else if ((reg [1] == "ridge") | (reg [1] == "lasso") | (reg [1] == "elastic"))
      {
        palpha = alpha
        if (reg [1] == "ridge")
          palpha = 0
        else if (reg [1] == "lasso")
          palpha = 1
        cv = NULL
        for (i in 1:10)
          cv = rbind (cv, glmnet::cv.glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda, standardize = TRUE, nfolds = 10)$cvm)
        cv = apply (cv, 2, mean)
        lambda.min = (lambda [101:1]) [which.min (cv)]
        if (graph)
        {
          plotmsep (lambda, cv)
          plotcoeff (lambda, glmnet::glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda, standardize = FALSE), cv)
        }
        model = list (model = glmnet::glmnet (as.matrix (x), y, alpha = palpha, lambda = lambda.min, standardize = FALSE), method = "regularisation")
        class (model) = "model"
      }
      else if ((reg [1] == "pcr") | (reg [1] == "plsr"))
      {
        model = list (model = lmfact (x = x, y = y, reg = reg [1], regeval = regeval, scale = scale, graph = graph, ...), method = "MRV")
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
    model = get (reg) (y~., data = x, ncomp = nbc, scale = scale)
    return (model)
  }

#' @keywords internal
lmselect <-
  function (x, y, regeval = c ("r2", "bic", "adjr2", "cp"), graph = TRUE, ...)
  {
    coeff = 1
    rss = leaps::regsubsets (y~., x, nvmax = ncol (x), method = "exhaustive")
    if (regeval [1] == "rsq")
      regeval = "r2"
    if (regeval [1] == "cp")
      regeval [1] = "Cp"
    if (graph)
      graphics::plot (rss, scale = regeval [1])
    if (regeval [1] == "r2")
      regeval = "rsq"
    if (regeval [1] == "Cp")
    {
      coeff = -1
      regeval = "cp"
    }
    if (regeval [1] == "bic")
    {
      coeff = -1
    }
    s = summary (rss)
    best = which.max (coeff * s [[regeval [1]]])
    var = colnames (s$which) [s$which [best, ]][-1]
    var = paste (var, collapse = "+")
    expr = paste ("lm (y~", var, ", x)", sep = "")
    model = eval (parse (text = expr))
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
#' @param params Object containing the parameters. If given, it replaces \code{size} and \code{decay}.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[nnet]{nnet}}
#' @examples
#' require (datasets)
#' data (trees)
#' MLPREG (trees [, -3], trees [, 3])
MLPREG <-
  function (x,
            y,
            size = 2:(ifelse (is.vector (x), 2, ncol (x))),
            decay = 10^(-3:-1),
            params = NULL,
            tune = FALSE,
            ...)
  {
    model = NULL
    d = cbind.data.frame (Y = y, x)
    minimum = apply (d, 2, min)
    range = apply (d, 2, max) - apply (d, 2, min)
    d.norm = sweep (sweep (d, 2, minimum, FUN = "-"), 2, range, FUN = "/")
    if (is.null (colnames (x)))
      colnames (d.norm) = c ("Y", paste ("X", 1:(ncol (d.norm) - 1), sep = ""))
    else
      colnames (d.norm) = c ("Y", colnames (x))
    if (!is.null (params))
    {
      size = params$hidden
      decay = params$decay
    }
    if (length (size) > 1 | length (decay) > 1)
    {
      tunecontrol = e1071::tune.control (sampling = "bootstrap", nboot = 20, boot.size = 1)
      model = e1071::tune.nnet (Y~., data = d.norm, size = size, decay = decay, tunecontrol = tunecontrol, ...)$best.model
    }
    else
      model = nnet::nnet (Y~., data = d.norm, size = size, decay = decay, trace = FALSE, ...)
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
  function (x, y, reg, regeval = c ("r2", "msep"), graph = TRUE, ...)
  {
    eval = toupper (regeval [1])
    optim = 0
    if (eval [1] == "R2")
      optim = 1
    else if (eval [1] == "MSEP")
      optim = -1
    res = utils::getFromNamespace (eval, ns = "pls") (utils::getFromNamespace (reg [1], ns = "pls") (y~., data = x, ncomp = min (ncol (x), nrow (x) - 2), validation = "LOO"), estimate = c ("train", "CV"))
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
      graphics::legend (pos, c ("Learning set", "LOOCV"), lty = 1:2, col = 1:2)
      graphics::abline (v = ncomp, lty = 2)
      bottom = min (res$val ["CV",,])
      top = max (res$val ["CV",,])
      delta = top - bottom
      graphics::text (ncomp, bottom + delta * .5, paste ("Nb. comp. :", ncomp), pos = 4)
    }
    return (ncomp)
  }

#' @keywords internal
plotcoeff <-
  function (lambda, model, cv = NULL)
  {
    x = log (lambda [101:1])
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
    x = log (lambda [101:1], base = 10)
    y = cv
    graphics::plot (x, y, t = "l", lwd = 2, xlab = expression(paste("log(", lambda, ")")), ylab = "MSEP",
          cex.lab = 1.5, cex.axis = 1.5, col = "red")
    bestx = x [which.min (y)]
    graphics::abline (v = bestx, lwd = 2, lty = 2, col = "darkgrey")
    graphics::text (bestx, max (y) * .95, bquote(paste("log(",lambda, ")=", .(bestx))), pos = 2, cex = 1.5)
    graphics::text (bestx, max (y) * .85, bquote(paste(lambda, "=", .(lambda [101:1] [which.min (y)]))), pos = 2, cex = 1.5)
  }

#' Polynomial Regression
#'
#' This function builds a polynomial regression model.
#' @name POLYREG
#' @param x Predictor \code{matrix}.
#' @param y Response \code{vector}.
#' @param degree The polynom degree.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{\link{model-class}}.
#' @export
#' @seealso \code{\link[mda]{polyreg}}
#' @examples
#' require (datasets)
#' data (trees)
#' POLYREG (trees [, -3], trees [, 3])
POLYREG <-
  function (x, y, degree = 2, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
      res = mda::polyreg (x, y, degree = degree)
    return (res)
  }

#' Plot function for a regression model
#'
#' Plot a regresion model on a 2-D plot. The predictor \code{x} should be one-dimensional.
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
    deltax = (max (x) - min (x)) * margin
    xlim = c (min (x) - deltax, max (x) + deltax)
    deltay = (max (y) - min (y)) * margin
    ylim = c (min (y) - deltay, max (y) + deltay)
    xl = data.frame (X = seq (xlim [1], xlim [2], length = 1000))
    if (!is.null (model$terms))
      colnames (xl) = attr (attr (model$terms, "factor"), "dimnames") [[1]] [2]
    graphics::plot (x, y, xaxs = "i", xlim = xlim, ylim = ylim, ...)
    graphics::lines (cbind (xl, stats::predict (model, xl)), col = 2)
  }

#' Plot the studentized residuals of a linear regression model
#'
#' Plot the studentized residuals of a linear regression model.
#'
#' @name resplot
#' @param model The model to be plotted.
#' @param index The index of the variable used for for the x-axis.
#' @export
#' @examples
#' require (datasets)
#' data (trees)
#' model = LINREG (trees [, -3], trees [, 3])
#' resplot (model) # Ordered by index
#' resplot (model, index = 0) # Ordered by variable "Volume" (dependant variable)
#' resplot (model, index = 1) # Ordered by variable "Girth" (independant variable)
#' resplot (model, index = 2) # Ordered by variable "Height" (independant variable)
resplot <-
  function (model, index = NULL)
  {
    mod = model
    xlab = "Index"
    if (methods::is (mod, "model"))
      mod = mod$model
    y = stats::rstudent (mod)
    if (is.null (index))
      index = 1:length (mod$residuals)
    else
    {
      xlab = colnames (mod$model) [index + 1]
      index = unlist (mod$model [index + 1])
    }
    ylim = c (min (y, -2) * 1.11, max (y, 2) * 1.11)
    graphics::plot (index, y, ylim = ylim, ylab = "Residuals", xlab = xlab, cex = 2, lwd = 2, col = ifelse (abs (y) <= 2, "darkgray", 2), cex.axis = 1.5, cex.lab = 1.5)
    mod = stats::loess (y ~ index)
    x = seq (min (index), max (index), length.out = 1001)
    graphics::lines (x, stats::predict (mod, x), lwd = 2, lty = 4, col = 4)
    graphics::abline (h = c (-2, 0, 2), lty = c (2, 1, 2), lwd = 2)
    select = which (abs (y) > 2)
    if (length (select) > 0)
      graphics::text (index [select], y [select], select, pos = 3, cex = 1)
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
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param params Object containing the parameters. If given, it replaces \code{epsilon}, \code{gamma} and \code{cost}.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVRl}}, \code{\link{SVRr}}
#' @examples
#' require (datasets)
#' data (trees)
#' SVR (trees [, -3], trees [, 3], kernel = "linear", cost = 1)
#' SVR (trees [, -3], trees [, 3], kernel = "radial", gamma = 1, cost = 1)
SVR <-
  function (x,
            y,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            kernel = c ("radial", "linear"),
            epsilon = c (.1, .5, 1),
            params = NULL,
            tune = FALSE,
            ...)
  {
    model = NULL
    if (!is.null (params))
    {
      gamma = params$gamma
      cost = params$cost
    }
    if (kernel [1] == "linear")
      gamma = 0
    if (length (gamma) > 1 | length (cost) > 1 | length (epsilon) > 1)
    {
      tunecontrol = e1071::tune.control(sampling = "bootstrap", nboot = 20, boot.size = 1)
      model = e1071::tune.svm (x, y, epsilon = epsilon,
                               gamma = gamma, cost = cost, kernel = kernel [1],
                               tunecontrol = tunecontrol, ...)$best.model
    }
    else
      model = e1071::svm (x, y, epsilon = c (.1, .5, 1), gamma = gamma, cost = cost, kernel = kernel [1], ...)
    res = model
    if (tune)
    {
      res = list (epsilon = model$epsilon, gamma = model$gamma, cost = model$cost)
      class (res) = "params"
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
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param params Object containing the parameters. If given, it replaces \code{epsilon}, \code{gamma} and \code{cost}.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVR}}
#' @examples
#' require (datasets)
#' data (trees)
#' SVRl (trees [, -3], trees [, 3], cost = 1)
SVRl <-
  function (x,
            y,
            cost = 2^(-3:3),
            epsilon = c (.1, .5, 1),
            params = NULL,
            tune = FALSE,
            ...)
  {
    return (SVR (
      x = x,
      y = y,
      cost = cost,
      epsilon = epsilon,
      kernel = "linear",
      params = params,
      tune = tune,
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
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param params Object containing the parameters. If given, it replaces \code{epsilon}, \code{gamma} and \code{cost}.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVR}}
#' @examples
#' require (datasets)
#' data (trees)
#' SVRr (trees [, -3], trees [, 3], gamma = 1, cost = 1)
SVRr <-
  function (x,
            y,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            epsilon = c (.1, .5, 1),
            params = NULL,
            tune = FALSE,
            ...)
  {
    return (SVR (
      x = x,
      y = y,
      gamma = gamma,
      cost = cost,
      epsilon = epsilon,
      kernel = "radial",
      params = params,
      tune = tune,
      ...
    ))
  }
