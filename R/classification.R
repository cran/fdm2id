#' Boosting methods model
#'
#' This class contains the classification model obtained by the CDA method.
#' @name boosting-class
#' @slot models List of models.
#' @slot x The learning set.
#' @slot y The target values.
#' @exportClass boosting
#' @seealso \code{\link{ADABOOST}}, \code{\link{BAGGING}}, \code{\link{predict.boosting}}
setClass ("boosting",
          representation (models = "list",
                          x = "data.frame",
                          y = "ANY"))

#' Canonical Disciminant Analysis model
#'
#' This class contains the classification model obtained by the CDA method.
#' @name cda-class
#' @slot proj The projection of the dataset into the canonical base. A \code{data.frame}.
#' @slot transform The transformation matrix between. A \code{matrix}.
#' @slot centers Coordinates of the class centers. A \code{matrix}.
#' @slot within The intr-class covarianc matrix. A \code{matrix}.
#' @slot eig The eigen-values. A \code{matrix}.
#' @slot dim The number of dimensions of the canonical base (numeric value).
#' @slot nb.classes The number of clusters (numeric value).
#' @slot train The training set (description). A \code{data.frame}.
#' @slot labels Class labels of the training set. Either a \code{factor} or an integer \code{vector}.
#' @slot model The prediction model.
#' @exportClass cda
#' @seealso \code{\link{CDA}}, \code{\link{plot.cda}}, \code{\link{predict.cda}}
setClass ("cda",
          representation (proj = "data.frame",
                          transform = "matrix",
                          centers = "matrix",
                          within = "matrix",
                          eig = "matrix",
                          dim = "numeric",
                          nb.classes = "numeric",
                          train = "data.frame",
                          labels = "factor",
                          model = "ANY"))

#' Training set and test set
#'
#' This class contains a dataset divided into four parts: the training set and test set, description and class labels.
#' @name dataset-class
#' @slot train.x the training set (description), as a \code{data.frame} or a \code{matrix}.
#' @slot train.y the training set (target), as a \code{vector} or a \code{factor}.
#' @slot test.x the training set (description), as a \code{data.frame} or a \code{matrix}.
#' @slot test.y the training set (target), as a \code{vector} or a \code{factor}.
#' @exportClass dataset
#' @seealso \code{\link{splitdata}}
setClass ("dataset",
          representation (train.x = "ANY",
                          train.y = "ANY",
                          test.x = "ANY",
                          test.y = "ANY"))

#' K Nearest Neighbours model
#'
#' This class contains the classification model obtained by the k-NN method.
#' @name knn-class
#' @slot train The training set (description). A \code{data.frame}.
#' @slot labels Class labels of the training set. Either a \code{factor} or an integer \code{vector}.
#' @slot k The \code{k} parameter.
#' @exportClass knn
#' @seealso \code{\link{KNN}}, \code{\link{predict.knn}}
setClass ("knn",
          representation (train = "data.frame",
                          labels = "factor",
                          k = "numeric"))

#' Generic classification or regression model
#'
#' This is a wrapper class containing the classification model obtained by any classification or regression method.
#' @name model-class
#' @slot model The wrapped model.
#' @slot method The name of the method.
#' @exportClass model
#' @seealso \code{\link{predict.model}}, \code{\link[stats]{predict}}
setClass ("model",
          representation (model = "ANY",
                          method = "character"))

#' Learning Parameters
#'
#' This class contains main parameters for various learning methods.
#' @name params-class
#' @slot decay The decay parameter.
#' @slot hidden The number of hidden nodes.
#' @slot epsilon The epsilon parameter.
#' @slot gamma The gamma parameter.
#' @slot cost The cost parameter.
#' @exportClass params
#' @seealso \code{\link{MLP}}, \code{\link{MLPREG}}, \code{\link{SVM}}, \code{\link{SVR}}
setClass ("params",
          representation (decay = "numeric",
                          hidden = "numeric",
                          epsilon = "numeric",
                          gamma = "numeric",
                          cost = "numeric"))

#' @keywords internal
adaboost.m1 <-
  function (x, y, learningmethod, nsamples, seed = NULL, ...)
  {
    set.seed (seed)
    if (is.vector (x))
      x = matrix (x, ncol = 1)
    w = rep (1 / nrow (x), nrow (x))
    epsilon = 0
    models = NULL
    iteration = 0
    while ((epsilon < .5) & (iteration < nsamples))
    {
      iteration = iteration + 1
      prob = w / sum (w)
      s = sample (nrow (x), nrow (x), replace = TRUE, prob = prob)
      xx = x [s, ]
      yy = y [s]
      model = learningmethod (xx, yy, ...)
      model$boostx = xx
      model$boosty = yy
      model$boostprob = prob
      rho = ifelse (predict (model, x) == y, 0, 1)
      epsilon = sum (prob * rho)
      if (epsilon > 0)
      {
        beta = epsilon / (1 - epsilon)
        model$boostweight = log (1 / beta)
        w = w * (beta^(1-rho))
        if (beta < 1)
          models = c (models, list (model))
      }
    }
    res = list (models = models, x = x, y = y)
    class (res) = "boosting"
    return (res)
  }

#' @keywords internal
adaboost.m2 <-
  function (x, y, learningmethod, nsamples, seed = NULL, ...)
  {
    set.seed (seed)
    if (is.vector (x))
      x = matrix (x, ncol = 1)
    Y = cbind (1:nrow (x), as.numeric (y))
    k = nlevels (y)
    D = rep (1, k * nrow (x))
    w = matrix (D / (k - 1), ncol = k)
    w [Y] = 0
    epsilon = 0
    models = NULL
    iteration = 0
    while ((epsilon < .5) & (iteration < nsamples))
    {
      iteration = iteration + 1
      W = apply (w, 1, sum)
      q = sweep (w, 1, W, "/")
      D = W / sum (W)
      s = sample (nrow (x), nrow (x), replace = TRUE, prob = D)
      xx = x [s, ]
      yy = y [s]
      model = learningmethod (xx, yy, ...)
      if (!is.null (model))
      {
        model$boostx = xx
        model$boosty = yy
        model$boostprob = D
        pred = predict (model, x, fuzzy = TRUE)
        if (!any (is.na (pred)))
        {
          rho = 1 - pred [Y] + apply (q * pred, 1, sum)
          epsilon = .5 * sum (D * rho)
          if (epsilon > 0)
          {
            beta = epsilon / (1 - epsilon)
            model$boostweight = log (1 / beta)
            exponent = .5 * (1 + sweep (-pred, 1, pred [Y], "+"))
            w = w * beta^exponent
            if (beta < 1)
              models = c (models, list (model))
          }
        }
      }
    }
    res = list (models = models, x = x, y = y)
    class (res) = "boosting"
    return (res)
  }

#' Classification using AdaBoost
#'
#' Ensemble learning, through AdaBoost Algorithm.
#' @name ADABOOST
#' @param x The dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' @param y The target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param learningmethod The boosted method.
#' @param nsamples The number of samplings.
#' @param fuzzy Indicates whether or not fuzzy classification should be used or not.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param seed A specified seed for random number generation.
#' @param ... Other specific parameters for the leaning method.
#' @return The classification model.
#' @export
#' @seealso \code{\link{BAGGING}}, \code{\link{predict.boosting}}
#' @examples
#' require (datasets)
#' data (iris)
#' ADABOOST (iris [, -5], iris [, 5], NB)
ADABOOST <-
  function (x, y, learningmethod, nsamples = 100, fuzzy = FALSE, tune = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.factor (y))
      {
        if (fuzzy)
          res = adaboost.m2 (x, y, learningmethod, nsamples, seed, ...)
        else
          res = adaboost.m1 (x, y, learningmethod, nsamples, seed, ...)
      }
      else
        res = NULL
    }
    return (res)
  }

#' Classification using Bagging
#'
#' Ensemble learning, through Bagging Algorithm.
#' @name BAGGING
#' @param x The dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' @param y The target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param learningmethod The boosted method.
#' @param nsamples The number of samplings.
#' @param size The size of the samples.
#' @param seed A specified seed for random number generation.
#' @param ... Other specific parameters for the leaning method.
#' @return The classification model.
#' @export
#' @seealso \code{\link{ADABOOST}}, \code{\link{predict.boosting}}
#' @examples
#' require (datasets)
#' data (iris)
#' BAGGING (iris [, -5], iris [, 5], NB)
BAGGING <-
  function (x, y, learningmethod, nsamples = 100, size = nrow (x), seed = NULL, ...)
  {
    set.seed (seed)
    if (is.vector (x))
      x = matrix (x, ncol = 1)
    s = matrix (sample (nrow (x), nsamples * size, replace = TRUE), ncol = nsamples)
    models = apply (s, 2, function (v)
    {
      train = x [v, ]
      target = y [v]
      model = learningmethod (train, target, ...)
      model$boostweight = 1
      return (model)
    })
    res = list (models = models, x = x, y = y)
    class (res) = "boosting"
    return (res)
  }

#' Bootstrap evaluation
#'
#' Evaluation a classification or regression method using bootstrap approach.
#' @name bootstrap
#' @param methods The classification or regression method to be evaluated.
#' @param x The dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' @param y The target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param eval The evaluation function.
#' @param nruns The number of bootstrap runs.
#' @param methodparameters Method parameters (if null tuning is done by cross-validation).
#' @param names Method names.
#' @param seed A specified seed for random number generation (useful for testing different method with the same bootstap samplings).
#' @param ... Other specific parameters for the leaning method.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluate}}, \code{\link{evaluation}}, \code{\link{bootstrap.curves}}
#' @examples
#' require ("datasets")
#' data (iris)
#' # One method, one evaluation criterion
#' bootstrap (NB, iris [, -5], iris [, 5], seed = 0)
#' # One method, two evaluation criteria
#' bootstrap (NB, iris [, -5], iris [, 5], eval = c ("accuracy", "kappa"), seed = 0)
#' # Three methods, two evaluation criteria
#' bootstrap (c (NB, LDA, LR), iris [, -5], iris [, 5], eval = c ("accuracy", "kappa"), seed = 0)
#' # List of methods in a variable
#' classif = c (NB, LDA, LR)
#' bootstrap (classif, iris [, -5], iris [, 5], eval = c ("accuracy", "kappa"), seed = 0,
#'            names = c ("NB", "LDA", "LR"))
#' # List of strings (method names)
#' classif = c ("NB", "LDA", "LR")
#' bootstrap (classif, iris [, -5], iris [, 5], eval = c ("accuracy", "kappa"), seed = 0)
bootstrap <-
  function (methods, x, y, eval = ifelse (is.factor (y), "accuracy", "r2"), nruns = 10, seed = NULL, methodparameters = NULL, names = NULL, ...)
  {
    set.seed (seed)
    methodNames = names
    if (is.character (methods))
    {
      methodNames = methods
      methods = sapply (methods, get)
    }
    else
    {
      if (is.null (names))
      {
        methodNames = as.character (match.call ()$methods)
        if (length (methodNames) > 1)
          methodNames = methodNames [-1]
        if (length (methodNames) != length (methods))
          methodNames = NULL
      }
    }
    if (is.vector (x))
      x = data.frame (X = x)
    predictions = NULL
    targets = NULL
    n = length (y)
    indices = 1:length (methods)
    if (is.null (methodparameters))
    {
      if (length (methods) == 1)
        methodparameters = methods (x, y, tune = TRUE, ...)
      else
        methodparameters = sapply (methods, function (method) method (x, y, tune = TRUE, ...))
    }
    set.seed (seed)
    samples = matrix (sample (n, n * nruns, replace = TRUE), ncol = nruns)
    for (i in 1:nruns)
    {
      s = samples [, i]
      targets = c (targets, y [-s])
      learn = x [s, ]
      if (is.vector (learn))
        learn = data.frame (X = learn)
      test = x [-s, ]
      if (is.vector (test))
        test = data.frame (X = test)
      rownames (learn) = 1:nrow (learn)
      if (length (methods) == 1)
      {
        models = methods (learn, y [s], graph = FALSE, methodparameters = methodparameters, ...)
        predictions = c (predictions, stats::predict (models, test, ...))
      }
      else
      {
        models = lapply (indices, function (i) methods [[i]] (learn, y [s], graph = FALSE, methodparameters = methodparameters [[i]], ...))
        predictions = rbind (predictions, sapply (models, function (model) as.numeric (stats::predict (model, test, ...))))
      }
    }
    if (is.factor (y))
    {
      lab = levels (y)
      l1 = length (lab)
      if (length (methods) == 1)
      {
        l2 = length (unique (predictions))
        if (l2 > l1)
          lab = c (lab, rep ("Unknown", l2 - l1))
        lab = lab [sort (unique (predictions))]
        predictions = factor (predictions, labels = lab)
      }
      else
      {
        predictions = as.data.frame (predictions)
        predictions = lapply (predictions, function (column) {
          l2 = length (unique (column))
          lab2 = lab
          if (l2 > l1)
            lab2 = c (lab, rep ("Unknown", l2 - l1))
          lab2 = lab2 [sort (unique (column))]
          return (factor (column, labels = lab2))
        })
        predictions = as.data.frame (predictions)
      }
      targets = factor (targets, labels = levels (y))
    }
    res = NULL
    if (length (methods) == 1)
      res = evaluation (predictions = predictions, targets = targets, eval = eval, ...)
    else
    {
      res = t (as.data.frame (lapply (predictions, function (column) evaluation (predictions = column, targets = targets, eval = eval, ...))))
      rownames (res) = methodNames
    }
    return (res)
  }

#' Plot evaluation curves with bootstrap sampling
#'
#' Evaluation a classification method according to ROC Curves or Cost Curves using bootstrap approach.
#' @name bootstrap.curves
#' @param methods The classification or regression method to be evaluated.
#' @param x The dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' @param y The target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param nruns The number of bootstrap runs.
#' @param seed A specified seed for random number generation (useful for testing different method with the same bootstap samplings).
#' @param curve A character string indicating the type of curve to be plotted.
#' @param methodparameters Method parameters (if null tuning is done by cross-validation).
#' @param new A logical value indicating whether a new plot should be be created or not.
#' @param lty The line type (and color) specified as an integer.
#' @param names Method names.
#' @param ... Other specific parameters for the leaning method.
#' @export
#' @seealso \code{\link{bootstrap}}, \code{\link[ROCR]{prediction}}, \code{\link[ROCR]{performance}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' # One method
#' bootstrap.curves (NB, d [, -5], d [, 5], seed = 0)
#' # Three methods
#' bootstrap.curves (c (NB, LDA, LR), d [, -5], d [, 5], seed = 0)
bootstrap.curves <-
  function (methods, x, y, nruns = 10, seed = NULL, curve = c ("ROC", "Cost"), methodparameters = NULL,
            new = TRUE, lty = 1, names = NULL, ...)
  {
    set.seed (seed)
    methodNames = names
    if (is.character (methods))
    {
      methodNames = methods
      methods = sapply (methods, get)
    }
    else
    {
      if (is.null (names))
      {
        methodNames = as.character (match.call ()$methods)
        if (length (methodNames) > 1)
          methodNames = methodNames [-1]
        if (length (methodNames) != length (methods))
          methodNames = NULL
      }
    }
    if (is.vector (x))
      x = data.frame (X = x)
    predictions = NULL
    targets = NULL
    n = length (y)
    indices = 1:length (methods)
    if (is.null (methodparameters))
    {
      if (length (methods) == 1)
        methodparameters = methods (x, y, tune = TRUE, ...)
      else
        methodparameters = sapply (methods, function (method) method (x, y, tune = TRUE, ...))
    }
    set.seed (seed)
    samples = matrix (sample (n, n * nruns, replace = TRUE), ncol = nruns)
    for (i in 1:nruns)
    {
      s = samples [, i]
      targets = c (targets, y [-s])
      learn = x [s, ]
      if (is.vector (learn))
        learn = data.frame (X = learn)
      test = x [-s, ]
      if (is.vector (test))
        test = data.frame (X = test)
      rownames (learn) = 1:nrow (learn)
      if (length (methods) == 1)
      {
        models = methods (learn, y [s], graph = FALSE, methodparameters = methodparameters [[i]], ...)
        predictions = c (predictions, stats::predict (models, test, ...))
      }
      else
      {
        models = lapply (indices, function (i) methods [[i]] (learn, y [s], graph = FALSE, methodparameters = methodparameters [[i]], ...))
        predictions = rbind (predictions, sapply (models, function (model) as.numeric (stats::predict (model, test, ...))))
      }
    }
    if (is.factor (y))
    {
      lab = levels (y)
      l1 = length (lab)
      if (length (methods) == 1)
      {
        l2 = length (unique (predictions))
        if (l2 > l1)
          lab = c (lab, rep ("Unknown", l2 - l1))
        lab = lab [sort (unique (predictions))]
        predictions = factor (predictions, labels = lab)
      }
      else
      {
        predictions = as.data.frame (predictions)
        predictions = lapply (predictions, function (column) {
          l2 = length (unique (column))
          lab2 = lab
          if (l2 > l1)
            lab2 = c (lab, rep ("Unknown", l2 - l1))
          lab2 = lab2 [sort (unique (column))]
          return (factor (column, labels = lab2))
        })
        predictions = as.data.frame (predictions)
      }
      targets = factor (targets, labels = levels (y))
    }
    res = NULL
    if (length (methods) == 1)
    {
      pred = ROCR::prediction (as.numeric (predictions), as.numeric (targets))
      type = tolower (curve [1])
      if (type == "roc")
      {
        perf = ROCR::performance (pred, "tpr", "fpr")
        ROCR::plot (perf, add = !new, lty = lty, col = lty, asp = 1)
      }
      if (type == "cost")
      {
        perf = ROCR::performance (pred, "ecost")
        ROCR::plot (perf, add = !new, lty = lty, col = lty)
      }
    }
    else
    {
      n = length (methods)
      add = c (!new, rep (T, n - 1))
      linetype = lty:(lty + n - 1)
      lapply (1:n, function (i) {
        pred = ROCR::prediction (as.numeric (predictions [, i]), as.numeric (targets))
        type = tolower (curve [1])
        if (type == "roc")
        {
          perf = ROCR::performance (pred, "tpr", "fpr")
          ROCR::plot (perf, add = add [i], lty = linetype [i], col = linetype [i], asp = 1)
        }
        if (type == "cost")
        {
          perf = ROCR::performance (pred, "ecost")
          ROCR::plot (perf, add = add [i], lty = linetype [i], col = linetype [i])
        }
        add = TRUE
        linetype = linetype + 1
      })
      graphics::legend ("bottomright", methodNames, lty = linetype, col = linetype)
    }
  }

#' Classification using CART
#'
#' This function builds a classification model using CART.
#' @name CART
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param minsplit The minimum leaf size during the learning.
#' @param maxdepth Set the maximum depth of any node of the final tree, with the root node counted as depth 0.
#' @param cp The complexity parameter of the tree. Cross-validation is used to determine optimal cp if NULL.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}, \code{\link{cartplot}}, \code{\link[rpart]{rpart}}
#' @examples
#' require (datasets)
#' data (iris)
#' CART (iris [, -5], iris [, 5])
CART <-
  function (train, labels, minsplit = 1, maxdepth = log2 (length (labels)), cp = NULL, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
      {
        train = matrix (train, ncol = 1)
        colnames (train) = "X"
      }
      d = cbind.data.frame (Class = labels, as.data.frame (train))
      complexity = cp
      if (is.null (cp))
      {
        model = rpart::rpart (Class~., d, minsplit = minsplit, xval = nrow (d), maxdepth = maxdepth, maxcompete = 0)
        mini = which.min (model$cptable [, 4])
        threshold = model$cptable [mini, 4] + model$cptable [mini, 5]
        complexity = model$cptable [which (model$cptable [, 4] < threshold) [1], 1]
      }
      model = rpart::rpart (Class~., d, minsplit = minsplit, cp = complexity, maxdepth = maxdepth, maxcompete = 0)
      res = list (model = model, method = "CART")
      class (res) = "model"
    }
    return (res)
  }

#' Depth
#'
#' Return the dept of a decision tree.
#' @name cartdepth
#' @param model The decision tree.
#' @return The depth.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}, \code{\link{cartplot}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartdepth (model)
cartdepth <-
  function (model)
  {
    return (ceiling (max (log (as.numeric (rownames (model$model$frame)), 2))) - 1)
  }

#' CART information
#'
#' Return various information on a CART model.
#' @name cartinfo
#' @param model The decision tree.
#' @return Various information organized into a vector.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartdepth}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}, \code{\link{cartplot}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartinfo (model)
cartinfo <-
  function (model)
  {
    return (c (Nodes = cartnodes (model), Leafs = cartleafs (model), Depth = cartdepth (model)))
  }

#' Number of Leafs
#'
#' Return the number of leafs of a decision tree.
#' @name cartleafs
#' @param model The decision tree.
#' @return The number of leafs.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartnodes}}, \code{\link{cartplot}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartleafs (model)
cartleafs <-
  function (model)
  {
    return (sum (model$model$frame$var == "<leaf>"))
  }

#' Number of Nodes
#'
#' Return the number of nodes of a decision tree.
#' @name cartnodes
#' @param model The decision tree.
#' @return The number of nodes.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartplot}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartnodes (model)
cartnodes <-
  function (model)
  {
    return (length (model$model$frame$var))
  }

#' CART Plot
#'
#' Plot a decision tree obtained by CART.
#' @name cartplot
#' @param model The decision tree.
#' @param margin an extra fraction of white space to leave around the borders of the tree. (Long labels sometimes get cut off by the default computation).
#' @param branch controls the shape of the branches from parent to child node. Any number from 0 to 1 is allowed. A value of 1 gives square shouldered branches, a value of 0 give V shaped branches, with other values being intermediate.
#' @param uniform if \code{TRUE}, uniform vertical spacing of the nodes is used; this may be less cluttered when fitting a large plot onto a page. The default is to use a non-uniform spacing proportional to the error in the fit.
#' @param fancy Logical. If \code{TRUE}, nodes are represented by ellipses (interior nodes) and rectangles (leaves) and labeled by yval. The edges connecting the nodes are labeled by left and right splits.
#' @param pretty an alternative to the minlength argument, see \code{\link[rpart]{labels.rpart}}.
#' @param fwidth Relates to option \code{fancy} and the width of the ellipses and rectangles. If \code{fwidth < 1} then it is a scaling factor (default = 0.8). If \code{fwidth > 1} then it represents the number of character widths (for current graphical device) to use.
#' @param fheight Relates to option \code{fancy} and the width of the ellipses and rectangles. If \code{fwidth < 1} then it is a scaling factor (default = 0.8). If \code{fwidth > 1} then it represents the number of character heights (for current graphical device) to use.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartplot (model)
cartplot <-
  function (model, margin = .2, branch = .3, uniform = TRUE, fancy = TRUE, pretty = TRUE, fwidth = 0, fheight = 0, ...)
  {
    graphics::plot (model$model, margin = margin, branch = branch, uniform = uniform, ...)
    graphics::text (model$model, fancy = fancy, pretty = pretty, fwidth = fwidth, fheight = fheight, ...)
  }

#' Classification using Canonical Discriminant Analysis
#'
#' This function builds a classification model using Canonical Discriminant Analysis.
#' @name CDA
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{glmnet}.
#' @export
#' @seealso \code{\link{plot.cda}}, \code{\link{predict.cda}}, \code{\link{cda-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' CDA (iris [, -5], iris [, 5])
CDA <-
  function (train, labels, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else if (length (unique (labels)) == nlevels (labels))
    {
      ll = factor (labels)
      m = scale (train, scale = FALSE)
      n = nrow (m)
      l = levels (ll)
      k = nlevels (ll)
      dim = min (k - 1, ncol (m))
      V = t (m) %*% m / n
      class = matrix (m [ll == l [1],], ncol = ncol (m))
      g = apply (class, 2, mean)
      B = nrow (class) * g %*% t (g)
      for (i in 2:k)
      {
        class = matrix (m [ll == l [i],], ncol = ncol (m))
        g = apply (class, 2, mean)
        B = B + (nrow (class) * g %*% t (g))
      }
      B = B / n
      s = eigen (solve (V) %*% B)
      t = s$vectors [, 1:dim]
      p = m %*% t
      W = V - B
      o = apply (matrix (m [ll == l [1], ], ncol = ncol (m)), 2, mean)
      for (i in 2:k)
        o = rbind (o, apply (matrix (m [ll == l [i], ], ncol = ncol (m)), 2, mean))
      colnames (o) = colnames (train)
      rownames (o) = l
      if (dim > 1)
      {
        colnames (p) = paste ("Can.", 1:dim)
        e = s$values [1:dim]
        e = cbind (e, 100 * e^2 / sum (e^2))
        e = cbind (e, cumsum (e [, 2]))
        colnames (t) = paste ("Can.", 1:dim)
        rownames (t) = colnames (train)
        colnames (e) = c ("eigenvalue", "percentage of variance", "cumulative percentage of variance")
        rownames (e) = paste ("Can.", 1:dim)
      }
      else
      {
        e = c (s$values [1], 100, 100)
        names (t) = colnames (train)
        names (e) = c ("eigenvalue", "percentage of variance", "cumulative percentage of variance")
      }
      prior = rep (1 / k, k)
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      model = MASS::lda (labels ~ ., as.data.frame (train), prior = prior)
      res = list (proj = as.data.frame (Re (p)),
                  transform = Re (t),
                  centers = Re (o),
                  eig = Re (e),
                  within = Re (W),
                  dim = dim,
                  nb.classes = k,
                  train = as.data.frame (train),
                  labels = ll,
                  model = model)
      class (res) = "cda"
    }
    else
      message (paste ("Missing classe(s):", levels (labels) [table (labels) == 0]))
    return (res)
  }

#' @keywords internal
cda.transform <-
  function (model, newdata) t (t (newdata) - apply (model$train, 2, mean)) %*% model$transform

#' Plot Cost Curves
#'
#' This function plots Cost Curves of several classification predictions.
#' @name cost.curves
#' @param methods.names The name of the compared methods (\code{vector}).
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param labels Actual labels of the dataset (\code{factor} or \code{vector}).
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{roc.curves}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' model.nb = NB (d [, -5], d [, 5])
#' model.lda = LDA (d [, -5], d [, 5])
#' pred.nb = predict (model.nb, d [, -5])
#' pred.lda = predict (model.lda, d [, -5])
#' cost.curves (c ("NB", "LDA"), cbind (pred.nb, pred.lda), d [, 5])
cost.curves <-
  function (methods.names, predictions, labels)
  {
    pred = ROCR::prediction (as.numeric (predictions [, 1]), as.numeric (labels))
    perf = ROCR::performance (pred, "ecost")
    ROCR::plot (perf)
    for (i in 2:ncol (predictions))
    {
      pred = ROCR::prediction (as.numeric (predictions [, i]), as.numeric (labels))
      perf = ROCR::performance (pred, "ecost")
      ROCR::plot (perf, add = TRUE, lty = i, col = i)
    }
    graphics::legend ("topleft", methods.names, lty = 1:ncol (predictions), col = 1:ncol (predictions))
  }

#' @keywords internal
emptyparams <-
  function ()
  {
    res = list ()
    class (res) = "params"
    return (res)
  }

#' @keywords internal
eval.accuracy <-
  function (predictions, targets, precision, recall, ...)
  {
    p = as.numeric (predictions)
    l = as.numeric (targets)
    sum (diag (table (p, l))) / length (l)
  }

#' @keywords internal
eval.fmeasure <-
  function (predictions, targets, precision = NULL, recall = NULL, beta = 1, positive = levels (targets) [1], ...)
  {
    if (is.null (precision))
      precision = evaluation.precision (predictions, targets, positive)
    if (is.null (recall))
      recall = evaluation.recall (predictions, targets, positive)
    res = (1 + beta * beta) * precision * recall / (beta * beta * precision + recall)
    return (res)
  }

#' @keywords internal
eval.fowlkesmallows <-
  function (predictions, targets, precision = NULL, recall = NULL, positive = levels (targets) [1], ...)
  {
    if (is.null (precision))
      precision = evaluation.precision (predictions, targets, positive)
    if (is.null (recall))
      recall = evaluation.recall (predictions, targets, positive)
    res = sqrt (precision * recall)
    return (res)
  }

#' @keywords internal
eval.goodness <-
  function (predictions, targets, beta = 1, precision = NULL, recall = NULL, positive = levels (targets) [1], ...)
  {
    if (is.null (precision))
      precision = evaluation.precision (predictions, targets, positive)
    if (is.null (recall))
      recall = evaluation.recall (predictions, targets, positive)
    res = (beta * precision + recall) / (beta + 1)
    return (res)
  }

#' @keywords internal
eval.jaccard <-
  function (predictions, targets, precision = NULL, recall = NULL, positive = levels (targets) [1], ...)
  {
    if (is.null (precision))
      precision = evaluation.precision (predictions, targets, positive)
    if (is.null (recall))
      recall = evaluation.recall (predictions, targets, positive)
    pr = precision * recall
    res = pr / (precision + recall - pr)
    return (res)
  }

#' @keywords internal
eval.kappa <-
  function (predictions, targets, precision = NULL, recall = NULL, ...)
  {
    p = as.numeric (predictions)
    l = as.numeric (targets)
    irr::kappa2 (cbind (p, l), weight = "equal")$value
  }

#' @keywords internal
eval.precision <-
  function (predictions, targets, precision, recall, positive = levels (targets) [1], ...)
  {
    return (precision)
  }

#' @keywords internal
eval.recall <-
  function (predictions, targets, precision, recall, positive = levels (targets) [1], ...)
  {
    return (recall)
  }


#' Evaluate several classication (or regression) methods
#'
#' Evaluation a classification or regression method using bootstrap approach.
#' @name evaluate
#' @param methods The classification or regression method to be evaluated.
#' @param dataset The dataset to be split (\code{data.frame} or \code{matrix}).
#' @param target The column index of the target variable (class label or response variable).
#' @param size The size of the training set (as an integer value).
#' @param names Method names.
#' @param eval The evaluation function.
#' @param seed A specified seed for random number generation.
#' @param ... Other specific parameters for the leaning method.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{bootstrap}}, \code{\link{evaluation}}, \code{\link{splitdata}}
#' @examples
#' require ("datasets")
#' data (iris)
#' evaluate (c (NB, LDA), iris, target = 5, eval = c ("accuracy", "kappa"), seed = 0)
evaluate <- function (methods, dataset, target = NULL, size = round (0.7 * nrow (dataset)), names = NULL, eval = "accuracy", seed = NULL, ...)
{
  methodNames = names
  if (is.character (methods))
  {
    methodNames = methods
    methods = sapply (methods, get)
  }
  else
  {
    if (is.null (names))
    {
      methodNames = as.character (match.call ()$methods)
      if (length (methodNames) > 1)
        methodNames = methodNames [-1]
      if (length (methodNames) != length (methods))
        methodNames = NULL
    }
  }
  d = NULL
  if ("dataset" %in% class (dataset))
    d = dataset
  else
  {
    if (!is.null (target))
      d = splitdata (dataset, target = target, size = size, seed = seed)
    else
      message ("Invalid dataset")
  }
  res = NULL
  if (length (methods) == 1)
    res = evaluation (predict (methods (d$train.x, d$train.y, ...), d$test.x, ...), d$test.y, eval = eval, ...)
  else
  {
    res = t (matrix (sapply (methods, function (method) evaluation (predict (method (d$train.x, d$train.y, ...), d$test.x, ...), d$test.y, eval = eval, ...)), ncol = length (methods)))
    colnames (res) = eval
    rownames (res) = methodNames
  }
  return (res)
}

#' Evaluation of classification or regression predictions
#'
#' Evaluation predictions of a classification or a regression model.
#' @name evaluation
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets The actual targets of the dataset (\code{factor} or \code{vector}).
#' @param eval The evaluation method.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation.msep}}, \code{\link{evaluation.r2}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' # Default evaluation for classification
#' evaluation (pred.nb, d$test.y)
#' # Evaluation with two criteria
#' evaluation (pred.nb, d$test.y, eval = c ("accuracy", "kappa"))
#' data (trees)
#' d = splitdata (trees, 3)
#' model.linreg = LINREG (d$train.x, d$train.y)
#' pred.linreg = predict (model.linreg, d$test.x)
#' # Default evaluation for regression
#' evaluation (pred.linreg, d$test.y)
evaluation <-
  function (predictions, targets, eval = ifelse (is.factor (targets), "accuracy", "r2"), ...)
  {
    precision = NULL
    recall = NULL
    if (is.factor (targets) & nlevels (targets) == 2)
    {
      precision = evaluation.precision (predictions = predictions, targets = targets, ...)
      recall = evaluation.recall (predictions, targets, ...)
    }
    res = NULL
    for (e in eval)
    {
      tmp = get (paste ("eval", e, sep = ".")) (predictions = predictions, targets = targets, precision = precision, recall = recall, ...)
      res = c (res, tmp)
    }
    if (!is.null (res))
      names (res) = eval
    return (res)
  }

#' Accuracy of classification predictions
#'
#' Evaluation predictions of a classification model according to accuracy.
#' @name evaluation.accuracy
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @return The evaluation of the predictions (numeric value).
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.accuracy (pred.nb, d$test.y)
evaluation.accuracy <-
  function (predictions, targets, ...)
  {
    return (eval.accuracy (predictions, targets, NULL, NULL))
  }

#' F-measure
#'
#' Evaluation predictions of a classification model according to the F-measure index.
#' @name evaluation.fmeasure
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param beta The weight given to precision.
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.fmeasure (pred.nb, d$test.y)
evaluation.fmeasure <-
  function (predictions, targets, beta = 1, positive = levels (targets) [1], ...)
  {
    return (eval.fmeasure (predictions, targets, beta = beta, positive = positive))
  }

#' Fowlkes–Mallows index
#'
#' Evaluation predictions of a classification model according to the Fowlkes–Mallows index.
#' @name evaluation.fowlkesmallows
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.fowlkesmallows (pred.nb, d$test.y)
evaluation.fowlkesmallows <-
  function (predictions, targets, positive = levels (targets) [1], ...)
  {
    return (eval.fowlkesmallows (predictions, targets, positive = positive))
  }

#' Goodness
#'
#' Evaluation predictions of a classification model according to Goodness index.
#' @name evaluation.goodness
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param beta The weight given to precision.
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.goodness (pred.nb, d$test.y)
evaluation.goodness <-
  function (predictions, targets, beta = 1, positive = levels (targets) [1], ...)
  {
    return (eval.goodness (predictions, targets, beta = beta, positive = positive))
  }

#' Jaccard index
#'
#' Evaluation predictions of a classification model according to Jaccard index.
#' @name evaluation.jaccard
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.jaccard (pred.nb, d$test.y)
evaluation.jaccard <-
  function (predictions, targets, positive = levels (targets) [1], ...)
  {
    return (eval.fmeasure (predictions, targets, positive = positive))
  }

#' Kappa evaluation of classification predictions
#'
#' Evaluation predictions of a classification model according to kappa.
#' @name evaluation.kappa
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}}, \code{\link{evaluation.precision}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.kappa (pred.nb, d$test.y)
evaluation.kappa <-
  function (predictions, targets, ...)
  {
    return (eval.kappa (predictions, targets))
  }

#' Precision of classification predictions
#'
#' Evaluation predictions of a classification model according to precision. Works only for two classes problems.
#' @name evaluation.precision
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}},
#' \code{\link{evaluation.recall}},\code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.precision (pred.nb, d$test.y)
evaluation.precision <-
  function (predictions, targets, positive = levels (targets) [1], ...)
  {
    if (nlevels (targets) != 2)
      stop ("evaluation.precision only works on data with two classes in the current implementation")
    t = table (targets, predictions)
    return (t [positive, positive] / sum (t [positive, ]))
  }

#' Recall of classification predictions
#'
#' Evaluation predictions of a classification model according to recall. Works only for two classes problems.
#' @name evaluation.recall
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param targets Actual targets of the dataset (\code{factor} or \code{vector}).
#' @param positive The label of the positive class.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' d = splitdata (d, 5)
#' model.nb = NB (d$train.x, d$train.y)
#' pred.nb = predict (model.nb, d$test.x)
#' evaluation.recall (pred.nb, d$test.y)
evaluation.recall <-
  function (predictions, targets, positive = levels (targets) [1], ...)
  {
    if (nlevels (targets) != 2)
      stop ("evaluation.precision only works on data with two classes in the current implementation")
    t = table (targets, predictions)
    return (t [positive, positive] / sum (t [, positive]))
  }

#' Classification using Gradient Boosting
#'
#' This function builds a classification model using Gradient Boosting
#' @name GRADIENTBOOSTING
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param ntree The number of trees in the forest.
#' @param learningrate The learning rate (between 0 and 1).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[xgboost]{xgboost}}
#' @examples
#' require (datasets)
#' data (iris)
#' GRADIENTBOOSTING (iris [, -5], iris [, 5])
GRADIENTBOOSTING <-
  function (train, labels,
            ntree = 500,
            learningrate = 0.3,
            tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      l = as.numeric (labels) - 1
      k = nlevels (labels)
      model = xgboost::xgboost (data = as.matrix (train), label = l, nrounds = ntree, objective = "multi:softprob", num_class = k, verbose = 0)
      res = list (model = model, lev = levels (labels), method = "XGB")
      class (res) = "model"
    }
    return (res)
  }

#' Classification using k-NN
#'
#' This function builds a classification model using Logistic Regression.
#' @name KNN
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param k The k parameter.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[class]{knn}}
#' @examples
#' require (datasets)
#' data (iris)
#' KNN (iris [, -5], iris [, 5])
KNN <-
  function (train, labels, k = 1:10, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      kk = k [1]
      if (is.vector (k) && (length (k) > 1))
      {
        tunecontrol = e1071::tune.control (sampling = "bootstrap", nboot = 20, boot.size = 1)
        kk = e1071::tune.knn (train, labels, k = k, tunecontrol = tunecontrol)$best.model$k
      }
      res = list (train = train, labels = labels, k = kk)
      class (res) = "knn"
    }
    return (res)
  }

#' Classification using Linear Discriminant Analysis
#'
#' This function builds a classification model using Linear Discriminant Analysis.
#' @name LDA
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[MASS]{lda}}
#' @examples
#' require (datasets)
#' data (iris)
#' LDA (iris [, -5], iris [, 5])
LDA <-
  function (train, labels, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      model = MASS::lda (x = train, grouping = labels)
      res = list (model = model, method = "LDA")
      class (res) = "model"
    }
    return (res)
  }

#' Classification using Logistic Regression
#'
#' This function builds a classification model using Logistic Regression.
#' @name LR
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[nnet]{multinom}}
#' @examples
#' require (datasets)
#' data (iris)
#' LR (iris [, -5], iris [, 5])
LR <-
  function (train, labels, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else if (length (unique (labels)) == nlevels (labels))
    {
      if (is.vector (train))
      {
        train = matrix (train, ncol = 1)
        colnames (train) = "X"
      }
      data = cbind.data.frame (train, Class = labels)
      model = nnet::multinom (formula = Class~., data, trace = FALSE)
      res = list (model = model, method = "LR")
      class (res) = "model"
    }
    else
      message (paste ("Missing classes:", levels (labels) [table (labels) == 0]))
    return (res)
  }

#' Classification using Multilayer Perceptron
#'
#' This function builds a classification model using Multilayer Perceptron.
#' @name MLP
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param size The size of the hidden layer (if a vector, cross-over validation is used to chose the best size).
#' @param decay The decay (between 0 and 1) of the backpropagation algorithm (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{size} and \code{decay}.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[nnet]{nnet}}
#' @examples
#' require (datasets)
#' data (iris)
#' MLP (iris [, -5], iris [, 5], size = 4, decay = .1)
MLP <-
  function (train,
            labels,
            size = ifelse (is.vector (train), 2:(1 + nlevels (labels)), 2:(ncol (train) + nlevels (labels))),
            decay = 10^(-3:-1),
            methodparameters = NULL,
            tune = FALSE,
            ...)
  {
    model = NULL
    if (is.vector (train))
      train = data.frame (X = train)
    d = cbind.data.frame (Class = labels, train)
    if (!is.null (methodparameters))
    {
      size = methodparameters$hidden
      decay = methodparameters$decay
    }
    if (length (size) > 1 | length (decay) > 1)
    {
      tunecontrol = e1071::tune.control(sampling = "bootstrap", nboot = 20, boot.size = 1)
      model = e1071::tune.nnet (Class~., data = d, size = size, decay = decay,
                                tunecontrol = tunecontrol, ...)$best.model
    }
    else
      model = nnet::nnet (Class~., data = d, size = size, decay = decay, trace = FALSE, ...)
    res = NULL
    if (tune)
    {
      res = list (decay = model$decay, hidden = model$n [2])
      class (res) = "params"
    }
    else
    {
      res = list (model = model, method = "MLP")
      class (res) = "model"
    }
    return (res)
  }

#' Classification using Naive Bayes
#'
#' This function builds a classification model using Naive Bayes.
#' @name NB
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{naiveBayes}}
#' @examples
#' require (datasets)
#' data (iris)
#' NB (iris [, -5], iris [, 5])
NB <-
  function (train, labels, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      model = e1071::naiveBayes (train, labels)
      res = list (model = model, method = "NB")
      class (res) = "model"
    }
    return (res)
  }

#' Plot function for cda-class
#'
#' Plot the learning set (and test set) on the canonical axes obtained by Canonical Discriminant Analysis (function \code{CDA}).
#' @name plot.cda
#' @param x The classification model (object of class \code{cda-class}).
#' @param newdata The test set (\code{matrix} or \code{data.frame}).
#' @param axes The canonical axes to be printed (numeric \code{vector}).
#' @param ... Other parameters.
#' @method plot cda
#' @export
#' @seealso \code{\link{CDA}}, \code{\link{predict.cda}}, \code{\link{cda-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CDA (iris [, -5], iris [, 5])
#' plot (model)
plot.cda <-
  function (x, newdata = NULL, axes = 1:2, ...)
  {
    if (x$dim > 1)
    {
      graphics::layout (rbind (1, 2), heights = c (1, 7))
      X = x$proj
      col = x$labels
      n = nrow (X)
      pch = rep (1, n)
      mar = graphics::par ()$mar
      mar [3] = 1
      opar = graphics::par (mar = c (0, 4, 0, 3))
      on.exit (graphics::par (opar))
      graphics::plot.new ()
      labels = c ("Training set", "Class centers")
      lpch = c (1, 19)
      if (!is.null (newdata))
      {
        X = rbind (X, cda.transform (x, newdata))
        col = c (col, predict.cda (x, newdata))
        m = nrow (newdata)
        pch = c (pch, rep (3, m))
        labels = c (labels, "Test set")
        lpch = c (lpch, 3)
      }
      graphics::legend ("topright", labels, pch = lpch, bty = "n")
      graphics::legend ("topleft", levels (x$labels), col = 2:(nlevels (x$labels) + 1), lty = 1, bty = "n")
      xlab = paste (colnames (X) [axes [1]], " (", round (x$eig [axes [1], 2], 2), " %)", sep = "")
      ylab = paste (colnames (X) [axes [2]], " (", round (x$eig [axes [2], 2], 2), " %)", sep = "")
      graphics::par (mar = mar)
      graphics::plot (X [, axes], col = unclass (col) + 1, pch = pch, xlab = xlab, ylab = ylab, asp = 1)
      for (i in 1:nlevels (x$labels)) graphics::points (t (apply(x$proj [x$labels == levels (x$labels) [i], axes], 2, mean)), col = i + 1, pch = 19)
    }
    else
    {
      graphics::layout (rbind (1, 2), heights = c (1, 7))
      X = Re (x$proj [, 1])
      col = x$labels
      n = length (X)
      pch = rep (1, n)
      mar = graphics::par ()$mar
      mar [3] = 1
      opar = graphics::par (mar = c (0, 4, 0, 3))
      on.exit (graphics::par (opar))
      graphics::plot.new()
      labels = c ("Training set", "Class centers")
      lpch = c (1, 19)
      if (!is.null (newdata))
      {
        X = c (X, Re (cda.transform (x, newdata)))
        col = c (col, stats::predict (x, newdata))
        m = nrow (newdata)
        pch = c (pch, rep (3, m))
        labels = c (labels, "Test set")
        lpch = c (lpch, 3)
      }
      graphics::legend ("topright", labels, pch = lpch, bty = "n")
      graphics::legend ("topleft", levels (x$labels), col = 2:(nlevels (x$labels) + 1), lty = 1, bty = "n")
      xlab = "Index"
      ylab = "Can. 1"
      graphics::par (mar = mar)
      X = cbind (1:n, X)
      xx = X [1:n, ]
      graphics::plot (X, col = unclass (col) + 1, pch = pch, xlab = xlab, ylab = ylab)
      for (i in 1:nlevels (x$labels)) graphics::points (t (apply(xx [x$labels == levels (x$labels) [i], ], 2, mean)), col = i + 1, pch = 19)
    }
    graphics::layout (1)
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained by a boosting method.
#' @name predict.boosting
#' @param object The classification model (of class \code{\link{boosting-class}}, created by \code{\link{ADABOOST}} or \code{\link{BAGGING}}).
#' @param test The test set (a \code{data.frame})
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @return A vector of predicted values (\code{factor}).
#' @param ... Other parameters.
#' @export
#' @method predict boosting
#' @seealso \code{\link{ADABOOST}}, \code{\link{BAGGING}}, \code{\link{boosting-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = BAGGING (d$train.x, d$train.y, NB)
#' predict (model, d$test.x)
#' model = ADABOOST (d$train.x, d$train.y, NB)
#' predict (model, d$test.x)
predict.boosting <- function (object, test, fuzzy = FALSE, ...)
{
  pred = lapply (object$models, function (model) predict (model, test, fuzzy, ...))
  weights = sapply (object$models, function (model) model$boostweight)
  res = NULL
  if (is.factor (object$y))
  {
    if (fuzzy)
      res = apply (simplify2array (pred), 1:2, mean)
    else
    {
      pred = sapply (pred, function (v) v)
      pred = apply (pred, 1, function (v) names (which.max (questionr::wtd.table (v, weights = weights))))
      labels = levels (object$y)
      pred = c (labels, pred)
      res = factor (pred, levels = labels) [-(1:(length (labels)))]
    }
  }
  else
  {
    pred = sapply (pred, function (v) v)
    res = rowMeans (pred)
  }
  return (res)
}

#' Model predictions
#'
#' This function predicts values based upon a model trained by \code{\link{CDA}}.
#' @name predict.cda
#' @param object The classification model (of class \code{\link{cda-class}}, created by \code{\link{CDA}}).
#' @param test The test set (a \code{data.frame})
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @return A vector of predicted values (\code{factor}).
#' @param ... Other parameters.
#' @export
#' @method predict cda
#' @seealso \code{\link{CDA}}, \code{\link{plot.cda}}, \code{\link{cda-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = CDA (d$train.x, d$train.y)
#' predict (model, d$test.x)
predict.cda <-
  function (object, test, fuzzy = FALSE, ...)
  {
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    res = NULL
    if (fuzzy)
      res = stats::predict (object$model, as.data.frame (test))$posterior
    else
      res = stats::predict (object$model, as.data.frame (test))$class
    return (res)
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained by \code{\link{KNN}}.
#' @name predict.knn
#' @param object The classification model (of class \code{\link[class]{knn}}).
#' @param test The test set (a \code{data.frame}).
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @param ... Other parameters.
#' @return A vector of predicted values (\code{factor}).
#' @export
#' @method predict knn
#' @seealso \code{\link{KNN}}, \code{\link{knn-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = KNN (d$train.x, d$train.y)
#' predict (model, d$test.x)
predict.knn <-
  function (object, test, fuzzy = FALSE, ...)
  {
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    res = NULL
    if (fuzzy)
      res = attr (caret::knn3Train (object$train, test, object$labels, object$k, prob = TRUE), "prob")
    else
      res = class::knn (object$train, test, object$labels, object$k)
    return (res)
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained by any classification or regression model.
#' @name predict.model
#' @param object The classification model (of class \code{\link{cda-class}}, created by \code{\link{CDA}}).
#' @param test The test set (a \code{data.frame}).
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @param ... Other parameters.
#' @return A vector of predicted values (\code{factor}).
#' @export
#' @method predict model
#' @seealso \code{\link{model-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = LDA (d$train.x, d$train.y)
#' predict (model, d$test.x)
predict.model <-
  function (object, test, fuzzy = FALSE, ...)
  {
    res = NULL
    if (object$method == "CART")
    {
      if (is.vector (test))
      {
        test = matrix (test, ncol = 1)
        colnames (test) = "X"
      }
      if (fuzzy)
        res = stats::predict (object$model, as.data.frame (test))
      else
        res = stats::predict (object$model, as.data.frame (test), type = "class")
    }
    else if (object$method == "LDA")
    {
      if (is.vector (test))
        test = matrix (test, ncol = 1)
      if (fuzzy)
      {
        res = stats::predict (object$model, test)$posterior
        k = length (object$model$lev)
        if (ncol (res) < k)
        {
          tmp = matrix (rep (0, k * nrow (test)), ncol = k)
          colnames (tmp) = object$model$lev
          tmp [, colnames (res)] = res
          res = tmp
        }
      }
      else
        res = stats::predict (object$model, test)$class
    }
    else if (object$method == "lm")
    {
      if (is.vector (test))
        test = data.frame (X = test)
      res = stats::predict (object$model, test)
    }
    else if (object$method == "LR")
    {
      if (is.vector (test))
        test = matrix (test, ncol = 1)
      if (ncol (test) == 1)
        colnames (test) = "X"
      if (fuzzy)
      {
        res = stats::predict (object$model, test, type = "probs")
        labels = object$model$lev
        if (length (labels) == 2)
        {
          res = cbind (res, 1 - res)
          colnames (res) = labels
        }
      }
      else
        res = stats::predict (object$model, test, type = "class")
    }
    else if (object$method == "MLP")
    {
      if (is.vector (test))
        test = data.frame (X = test)
      if (fuzzy)
      {
        res = stats::predict (object$model, test, type = "raw")
        labels = object$model$lev
        if (length (labels) == 2)
        {
          res = cbind (res, 1 - res)
          colnames (res) = labels
        }
      }
      else
      {
        pred = stats::predict (object$model, test, type = "class")
        labels = object$model$lev
        pred = c (labels, pred)
        res = factor (pred, levels = labels) [-(1:(length (labels)))]
      }
    }
    else if (object$method == "MLPREG")
    {
      d.norm = sweep (sweep (as.data.frame (test), 2, object$model$minimum [-1], FUN = "-"), 2, object$model$range [-1], FUN = "/")
      colnames (d.norm) = attr (attr (object$model$model$terms, "factor"), "dimnames") [[1]] [-1]
      res = (stats::predict (object$model$model, d.norm) * object$model$range [1]) + object$model$minimum [1]
    }
    else if (object$method == "MRV")
    {
      res = stats::predict (object$model, test, ncomp = object$model$ncomp)
    }
    else if (object$method == "NB")
    {
      type = "class"
      if (fuzzy)
        type = "raw"
      if (is.vector (test))
        test = matrix (test, ncol = 1)
      res = stats::predict (object$model, test, type = type)
    }
    else if (object$method == "regularisation")
    {
      res = stats::predict (object$model, as.matrix (test))
    }
    else if (object$method == "RF")
    {
      if (is.vector (test))
        test = matrix (test, ncol = 1)
      res = stats::predict (object$model, test)
    }
    else if (object$method == "SVM")
    {
      if (fuzzy)
        res = attr (stats::predict (object$model, test, probability = TRUE), "probabilities")
      else
        res = stats::predict (object$model, test, probability = FALSE)
    }
    else if (object$method == "XGB")
    {
      res = predict (object$model, as.matrix (test), reshape = T)
      if (fuzzy)
        colnames (res) = object$lev
      else
      {
        levels = 1:length (object$lev)
        res = factor (c (levels, apply (res, 1, which.max)), labels = object$lev) [-levels]
      }
    }
    else
    {
      res = stats::predict (object$model, test, ...)
    }
    return (res)
  }

#' Classification using Random Forest
#'
#' This function builds a classification model using Random Forest
#' @name RANDOMFOREST
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param ntree The number of trees in the forest.
#' @param nvar Number of variables randomly sampled as candidates at each split.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[randomForest]{randomForest}}
#' @examples
#' require (datasets)
#' data (iris)
#' RANDOMFOREST (iris [, -5], iris [, 5])
RANDOMFOREST <-
  function (train, labels,
            ntree = 500,
            nvar = if (!is.null (labels) && !is.factor (labels)) max (floor (ncol (train)/3), 1) else floor (sqrt (ncol (train))),
            tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      res = randomForest::randomForest(x = train, y = labels, ntree = ntree, ntry = nvar, ...)
      res = list (model = res, method = "RF")
      class (res) = "model"
    }
    return (res)
  }

#' Plot ROC Curves
#'
#' This function plots ROC Curves of several classification predictions.
#' @name roc.curves
#' @param methods.names The name of the compared methods (\code{vector}).
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param labels Actual labels of the dataset (\code{factor} or \code{vector}).
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{cost.curves}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' model.nb = NB (d [, -5], d [, 5])
#' model.lda = LDA (d [, -5], d [, 5])
#' pred.nb = predict (model.nb, d [, -5])
#' pred.lda = predict (model.lda, d [, -5])
#' roc.curves (c ("NB", "LDA"), cbind (pred.nb, pred.lda), d [, 5])
roc.curves <-
  function (methods.names, predictions, labels)
  {
    pred = ROCR::prediction (as.numeric (predictions [, 1]), as.numeric (labels))
    perf = ROCR::performance (pred, "tpr", "fpr")
    ROCR::plot (perf)
    nclasses = ncol (predictions)
    for (i in 2:nclasses)
    {
      pred = ROCR::prediction (as.numeric (predictions [, i]), as.numeric (labels))
      perf = ROCR::performance (pred, "tpr", "fpr")
      ROCR::plot (perf, add = TRUE, lty = i, col = i)
    }
    graphics::legend ("bottomright", methods.names, lty = 1:nclasses, col = 1:nclasses)
  }

#' Classification using one-level decision tree
#'
#' This function builds a classification model using CART with maxdepth = 1.
#' @name STUMP
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param minsplit The minimum leaf size during the learning.
#' @param cp The complexity parameter of the tree. Cross-validation is used to determine optimal cp if NULL.
#' @param randomvar If true, the model uses a random variable.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link{CART}}
#' @examples
#' require (datasets)
#' data (iris)
#' STUMP (iris [, -5], iris [, 5])
STUMP <-
  function (train, labels, minsplit = 1, cp = NULL, randomvar = TRUE, tune = FALSE, ...)
  {
    new = train
    if (randomvar && (!is.vector (train)))
    {
      var = sample (ncol (train), 1)
      new = matrix (train [, var], ncol = 1)
      colnames (new) = colnames (train) [var]
    }
    return (CART (new, labels, minsplit = minsplit, 1, cp = cp, tune = tune, ...))
  }

#' Classification using Support Vector Machine
#'
#' This function builds a classification model using Support Vector Machine.
#' @name SVM
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param gamma The gamma parameter (if a vector, cross-over validation is used to chose the best size).
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param kernel The kernel type.
#' @param methodparameters Object containing the parameters. If given, it replaces \code{gamma} and \code{cost}.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVMl}}, \code{\link{SVMr}}
#' @examples
#' require (datasets)
#' data (iris)
#' SVM (iris [, -5], iris [, 5], kernel = "linear", cost = 1)
#' SVM (iris [, -5], iris [, 5], kernel = "radial", gamma = 1, cost = 1)
SVM <-
  function (train,
            labels,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            kernel = c ("radial", "linear"),
            methodparameters = NULL,
            tune = FALSE,
            ...)
  {
    model = NULL
    if (!is.null (methodparameters))
    {
      gamma = methodparameters$gamma
      cost = methodparameters$cost
    }
    if (kernel [1] == "linear")
      gamma = 0
    if (length (gamma) > 1 | length (cost) > 1)
    {
      tunecontrol = e1071::tune.control(sampling = "bootstrap", nboot = 20, boot.size = 1)
      model = e1071::tune.svm (train, labels,
                               gamma = gamma, cost = cost, kernel = kernel [1],
                               tunecontrol = tunecontrol, probability = TRUE, ...)$best.model
    }
    else
      model = e1071::svm (train, labels, gamma = gamma, cost = cost, kernel = kernel [1], probability = TRUE, ...)
    res = NULL
    if (tune)
    {
      res = list (gamma = model$gamma, cost = model$cost)
      class (res) = "params"
    }
    else
    {
      res = list (model = model, method = "SVM")
      class (res) = "model"
    }
    return (res)
  }

#' Classification using Support Vector Machine with a linear kernel
#'
#' This function builds a classification model using Support Vector Machine with a linear kernel.
#' @name SVMl
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{gamma} and \code{cost}.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVM}}
#' @examples
#' require (datasets)
#' data (iris)
#' SVMl (iris [, -5], iris [, 5], cost = 1)
SVMl <-
  function (train,
            labels,
            cost = 2^(-3:3),
            methodparameters = NULL,
            tune = FALSE,
            ...)
  {
    return (SVM (
      train = train,
      labels = labels,
      gamma = NULL,
      cost = cost,
      kernel = "linear",
      methodparameters = methodparameters,
      tune = tune,
      ...
    ))
  }

#' Classification using Support Vector Machine with a radial kernel
#'
#' This function builds a classification model using Support Vector Machine with a radial kernel.
#' @name SVMr
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param gamma The gamma parameter (if a vector, cross-over validation is used to chose the best size).
#' @param cost The cost parameter (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{gamma} and \code{cost}.
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVM}}
#' @examples
#' require (datasets)
#' data (iris)
#' SVMr (iris [, -5], iris [, 5], gamma = 1, cost = 1)
SVMr <-
  function (train,
            labels,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            methodparameters = NULL,
            tune = FALSE,
            ...)
  {
    return (SVM (
      train = train,
      labels = labels,
      gamma = gamma,
      cost = cost,
      kernel = "radial",
      methodparameters = methodparameters,
      tune = tune,
      ...
    ))
  }
