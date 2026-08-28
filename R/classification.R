#' Boosting methods model
#'
#' This class contains the ensemble of models obtained by a boosting or bagging method
#' (\code{\link{ADABOOST}}, \code{\link{BAGGING}}).
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{models}}{List of models.}
#'   \item{\code{x}}{The learning set.}
#'   \item{\code{y}}{The target values.}
#'   \item{\code{nsamples}}{The number of models that were asked for. Boosting keeps fewer
#'     when it runs out of models better than chance; \code{\link{print.boosting}} says so.}
#' }
#' @name boosting-class
#' @seealso \code{\link{ADABOOST}}, \code{\link{BAGGING}}, \code{\link{predict.boosting}}
NULL

#' Canonical Disciminant Analysis model
#'
#' This class contains the classification model obtained by the CDA method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{proj}}{The projection of the dataset into the canonical base. A \code{data.frame}.}
#'   \item{\code{transform}}{The transformation matrix between. A \code{matrix}.}
#'   \item{\code{centers}}{Coordinates of the class centers. A \code{matrix}.}
#'   \item{\code{within}}{The intra-class covariance matrix. A \code{matrix}.}
#'   \item{\code{eig}}{One row per canonical axis and four columns: the \code{eigenvalue}
#'     (of \eqn{V^{-1}B}, i.e. the squared canonical correlation, between 0 and 1), its
#'     \code{percentage of variance} (share of the trace) and the \code{cumulative} one, and
#'     the \code{discriminant power} -- the share of the trace of \eqn{W^{-1}B}, which is what
#'     \code{\link[MASS]{lda}} and most other software call the proportion of trace. A
#'     \code{matrix}, or a named \code{vector} when there is a single axis.}
#'   \item{\code{dim}}{The number of dimensions of the canonical base (numeric value).}
#'   \item{\code{nb.classes}}{The number of clusters (numeric value).}
#'   \item{\code{train}}{The training set (description). A \code{data.frame}.}
#'   \item{\code{labels}}{Class labels of the training set. Either a \code{factor} or an integer \code{vector}.}
#'   \item{\code{model}}{The prediction model.}
#' }
#' @name cda-class
#' @seealso \code{\link{CDA}}, \code{\link{plot.cda}}, \code{\link{predict.cda}}
NULL

#' Training set and test set
#'
#' This class contains a dataset divided into four parts: the training set and test set, description and class labels.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{train.x}}{the training set (description), as a \code{data.frame} or a \code{matrix}.}
#'   \item{\code{train.y}}{the training set (target), as a \code{vector} or a \code{factor}.}
#'   \item{\code{test.x}}{the training set (description), as a \code{data.frame} or a \code{matrix}.}
#'   \item{\code{test.y}}{the training set (target), as a \code{vector} or a \code{factor}.}
#' }
#' @name dataset-class
#' @seealso \code{\link{splitdata}}
NULL

#' K Nearest Neighbours model
#'
#' This class contains the classification model obtained by the k-NN method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{train}}{The training set (description). A \code{data.frame}.}
#'   \item{\code{labels}}{Class labels of the training set. Either a \code{factor} or an integer \code{vector}.}
#'   \item{\code{k}}{The \code{k} parameter.}
#' }
#' @name knn-class
#' @seealso \code{\link{KNN}}, \code{\link{predict.knn}}
NULL

#' Generic classification or regression model
#'
#' This is a wrapper class containing the classification model obtained by any classification or regression method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{model}}{The wrapped model.}
#'   \item{\code{method}}{The name of the method.}
#' }
#' @name model-class
#' @seealso \code{\link{predict.model}}, \code{\link[stats]{predict}}
NULL

#' Learning Parameters
#'
#' This class contains main parameters for various learning methods.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{decay}}{The decay parameter.}
#'   \item{\code{hidden}}{The number of hidden nodes.}
#'   \item{\code{epsilon}}{The epsilon parameter.}
#'   \item{\code{gamma}}{The gamma parameter.}
#'   \item{\code{cost}}{The cost parameter.}
#' }
#' @name params-class
#' @seealso \code{\link{MLP}}, \code{\link{MLPREG}}, \code{\link{SVM}}, \code{\link{SVR}}
NULL

#' @keywords internal
adaboost.m1 <-
  function (x, y, learningmethod, nsamples, seed = NULL, ...)
  {
    setseed (seed)
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
      # Not 'epsilon == 0': a model that is perfect up to rounding (LDA posteriors of
      # 1 - 1e-16 on well-separated classes) gives a minuscule but non-zero epsilon, hence a
      # beta so small that the next round of weights underflows to exactly zero -- and the
      # sampling probabilities become NaN.
      if (epsilon <= .Machine$double.eps^0.5)
      {
        # A base learner that classifies the whole reweighted training set correctly leaves
        # nothing to boost: it is kept, and the boosting stops there.
        model$boostweight = 1
        models = c (models, list (model))
        break
      }
      if (epsilon > 0)
      {
        beta = epsilon / (1 - epsilon)
        model$boostweight = log (1 / beta)
        w = w * (beta^(1-rho))
        if (beta < 1)
          models = c (models, list (model))
      }
    }
    return (boosting.result (models, x, y, nsamples, "ADABOOST"))
  }

#' @keywords internal
# Wraps the models an ensemble method ended up with, and says so when they are fewer than the
# caller asked for. Boosting stops as soon as a model is no better than chance, which is
# correct but was silent: ADABOOST (nsamples = 100) commonly kept four models on easy data
# with nothing to say so.
boosting.result <-
  function (models, x, y, nsamples, fname)
  {
    if (length (models) == 0)
      stop (fname, ": not a single model could be kept -- the very first one was no better ",
            "than chance on the training set, so there is nothing to boost. Try another base ",
            "method (a weak but non-trivial one, e.g. STUMP), or check the training data.")
    # Boosting ends as soon as a model is no better than chance on the reweighted training set,
    # or as soon as one is perfect on it, so an ensemble is often smaller than asked for. That
    # is normal and used to be announced by a message on every call; print.boosting() says it
    # instead, to whoever looks at the model.
    res = list (models = models, x = x, y = y, nsamples = nsamples)
    class (res) = "boosting"
    return (res)
  }

#' @keywords internal
adaboost.m2 <-
  function (x, y, learningmethod, nsamples, seed = NULL, ...)
  {
    setseed (seed)
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
      # Every weight of a row can underflow to zero after enough rounds, which turns q and the
      # sampling probabilities into NaN and stops sample() with "NA in probability vector".
      # There is nothing left to reweight at that point.
      if (any (!is.finite (W)) || (sum (W) == 0))
        break
      q = sweep (w, 1, W, "/")
      q [!is.finite (q)] = 0
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
          if (epsilon <= .Machine$double.eps^0.5)
          {
            # See adaboost.m1(): nothing left to boost, so the model is kept and the loop ends.
            model$boostweight = 1
            models = c (models, list (model))
            break
          }
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
    return (boosting.result (models, x, y, nsamples, "ADABOOST"))
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
#' @inheritParams tune.doc
#' @param seed A specified seed for random number generation.
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{ADABOOST} does not
#' support reusing pre-tuned parameters (the base learner given as \code{learningmethod} is
#' tuned independently on each boosting sample).
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{ADABOOST} does not produce a plot.
#' @param ... Other specific parameters for the leaning method.
#' @return The classification model.
#' @export
#' @seealso \code{\link{BAGGING}}, \code{\link{predict.boosting}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' ADABOOST (iris [, -5], iris [, 5], NB)
#' }
ADABOOST <-
  function (x, y, learningmethod, nsamples = 100, fuzzy = FALSE,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (!is.factor (y))
        stop ("ADABOOST: the target must be a factor -- AdaBoost is a classification method ",
              "and this one is numeric. Use BAGGING(), which handles both, or one of the ",
              "regression methods (GBREG, LINREG, SVR, ...).")
      if (fuzzy)
        res = adaboost.m2 (x, y, learningmethod, nsamples, seed, ...)
      else
        res = adaboost.m1 (x, y, learningmethod, nsamples, seed, ...)
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
#' @param bag.size The size of the samples.
#' @param seed A specified seed for random number generation.
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{BAGGING} does not
#' support reusing pre-tuned parameters (the base learner is fitted afresh on each sample).
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{BAGGING} does not produce a plot.
#' @param ... Other specific parameters for the leaning method.
#' @return The classification model.
#' @export
#' @seealso \code{\link{ADABOOST}}, \code{\link{predict.boosting}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' BAGGING (iris [, -5], iris [, 5], NB)
#' }
BAGGING <-
  function (x, y, learningmethod, nsamples = 100, bag.size = nrow (x),
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    # Explicit (and, for the last two, unused) formals, like every other learning function:
    # BAGGING forwards '...' to the base learner, so the arguments performance() always passes
    # have to be caught here rather than reach it.
    if (tune)
      return (emptyparams ())
    setseed (seed)
    if (is.vector (x))
      x = matrix (x, ncol = 1)
    s = matrix (sample (nrow (x), nsamples * bag.size, replace = TRUE), ncol = nsamples)
    models = apply (s, 2, function (v)
    {
      train = x [v, ]
      target = y [v]
      model = learningmethod (train, target, ...)
      model$boostweight = 1
      model$boostx = train
      model$boosty = target
      return (model)
    })
    return (boosting.result (models, x, y, nsamples, "BAGGING"))
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
#' @param xval The number of cross-validation folds used to choose \code{cp}, when \code{cp} is
#' \code{NULL}. \code{xval = nrow (train)} gives a leave-one-out cross-validation, which fits
#' one tree per observation and costs about \code{nrow (train) / 10} times as much.
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{CART} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{CART} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}, \code{\link{cartplot}}, \code{\link[rpart]{rpart}}
#' @examples
#' require (datasets)
#' data (iris)
#' CART (iris [, -5], iris [, 5])
CART <-
  function (train, labels, minsplit = 1, maxdepth = log2 (length (labels)), cp = NULL, xval = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
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
        # cp is chosen by cross-validation on 'xval' folds, then the tree refitted with it.
        model = rpart::rpart (Class~., d, minsplit = minsplit, xval = xval, maxdepth = maxdepth, maxcompete = 0, model = TRUE)
        mini = which.min (model$cptable [, 4])
        threshold = model$cptable [mini, 4] + model$cptable [mini, 5]
        complexity = model$cptable [which (model$cptable [, 4] < threshold) [1], 1]
      }
      model = rpart::rpart (Class~., d, minsplit = minsplit, cp = complexity, maxdepth = maxdepth, maxcompete = 0, model = TRUE)
      type = ifelse (is.factor (labels), "class", "reg")
      res = list (model = model, method = "CART", type = type)
      class (res) = "model"
    }
    return (res)
  }

#' Depth
#'
#' Return the depth of a decision tree.
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
    # rpart numbers its nodes so that the children of node i are 2i and 2i + 1: the depth of
    # a node is therefore floor (log2 (node)), and the depth of the tree is that of its
    # deepest node. The previous formula, ceiling (log2 (max)) - 1, happened to agree
    # everywhere except on a tree reduced to its root (max node = 1), where it returned -1.
    return (floor (log2 (max (as.numeric (rownames (model$model$frame))))))
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
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{CART}}, \code{\link{cartdepth}}, \code{\link{cartinfo}}, \code{\link{cartleafs}}, \code{\link{cartnodes}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CART (iris [, -5], iris [, 5])
#' cartplot (model)
cartplot <-
  function (model, ...)
  {
    col = as.list (sort (unique (as.numeric (model$model$model$Class) + 1)))
    rpart.plot::rpart.plot (model$model, box.palette = col, type = 0)
  }

#' Classification using Canonical Discriminant Analysis
#'
#' This function builds a classification model using Canonical Discriminant Analysis.
#'
#' The projection is computed from the class sizes, as the between-class scatter requires. The
#' predictions, on the other hand, use \emph{equal} prior probabilities -- an observation goes
#' to the nearest class centre in the canonical space, whatever the size of that class. This is
#' the geometric reading \code{\link{plot.cda}} draws, and it is where \code{CDA} differs from
#' \code{\link{LDA}}, which weights the classes by their observed frequencies: on an
#' imbalanced problem the two do not predict the same thing.
#' @name CDA
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{CDA} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{CDA} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{cda}.
#' @export
#' @seealso \code{\link{plot.cda}}, \code{\link{predict.cda}}, \code{\link{cda-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' CDA (iris [, -5], iris [, 5])
CDA <-
  function (train, labels, tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      # factor(labels) alone was not enough: the previous guard was
      # 'length (unique (labels)) == nlevels (labels)', which is never TRUE for a character
      # vector (nlevels() returns 0 there), so CDA silently returned NULL on perfectly valid
      # character labels. check.classes() handles factors and character vectors alike.
      ll = check.classes (labels, "CDA")
      labels = ll
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
      # Named 'tf', not 't': a local 't' shadows base::t() for the rest of the body, and the
      # transpositions above only kept working because they are evaluated before the
      # assignment.
      tf = s$vectors [, 1:dim, drop = FALSE]
      p = m %*% tf
      W = V - B
      o = apply (matrix (m [ll == l [1], ], ncol = ncol (m)), 2, mean)
      for (i in 2:k)
        o = rbind (o, apply (matrix (m [ll == l [i], ], ncol = ncol (m)), 2, mean))
      colnames (o) = colnames (train)
      rownames (o) = l
      eignames = c ("eigenvalue", "percentage of variance",
                    "cumulative percentage of variance", "discriminant power")
      if (dim > 1)
      {
        colnames (p) = paste ("Can.", 1:dim)
        e = Re (s$values [1:dim])
        # The share of the trace, computed on the very eigenvalues shown in the first column.
        e = cbind (e, 100 * e / sum (e))
        e = cbind (e, cumsum (e [, 2]))
        e = cbind (e, 100 * cda.power (Re (s$values [1:dim])))
        colnames (tf) = paste ("Can.", 1:dim)
        rownames (tf) = colnames (train)
        colnames (e) = eignames
        rownames (e) = paste ("Can.", 1:dim)
      }
      else
      {
        e = c (Re (s$values [1]), 100, 100, 100)
        # A single axis got no name at all, so as.data.frame() called it "V1" -- while
        # plot.cda() has always labelled it "Can. 1".
        colnames (p) = "Can. 1"
        tf = drop (tf)
        names (tf) = colnames (train)
        names (e) = eignames
      }
      prior = rep (1 / k, k)
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      model = MASS::lda (labels ~ ., as.data.frame (train), prior = prior)
      res = list (proj = as.data.frame (Re (p)),
                  transform = Re (tf),
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
    return (res)
  }

#' @keywords internal
cda.transform <-
  function (model, newdata) t (t (newdata) - apply (model$train, 2, mean)) %*% model$transform

#' @keywords internal
# The share of discriminant power carried by each canonical axis.
#
# CDA() diagonalises V^-1 B (V being the *total* covariance), so its eigenvalues are the
# squared canonical correlations, between 0 and 1: their share of the trace is what the
# "percentage of variance" column reports. Discriminant analysis is more usually read on the
# eigenvalues of W^-1 B (W being the *within-class* covariance), which measure between-class
# over within-class variance and are unbounded. The two are the same axes and are related by
#
#   lambda_W = lambda_V / (1 - lambda_V)
#
# so no second diagonalisation is needed. The share of trace computed on lambda_W is what
# MASS::lda(), SPSS and SAS call the proportion of trace: 99.1 / 0.9 on iris, where the
# "percentage of variance" column reads 81.4 / 18.6. Both describe the same two axes; the
# first says how much of the total dispersion they separate, the second how well they separate
# it.
#
# An axis that separates the classes perfectly has lambda_V = 1 and infinite discriminant
# power; such axes then share the whole of it between them.
cda.power <-
  function (values)
  {
    values = pmin (pmax (values, 0), 1)
    perfect = values >= 1
    if (any (perfect))
      return (perfect / sum (perfect))
    power = values / (1 - values)
    total = sum (power)
    if ((total == 0) || (!is.finite (total)))
      return (rep (1 / length (power), length (power)))
    return (power / total)
  }

#' Confusion matrix
#'
#' Plot a confusion matrix. Rows are the true labels and columns the predicted ones.
#' @name confusion
#' @param predictions The prediction. \strong{This is the first argument}, as for every other
#' evaluation function of the package (\code{\link{evaluation}}, \code{\link{evaluation.accuracy}},
#' \code{\link{evaluation.precision}}, ...): passing the ground truth first returns the transposed
#' matrix.
#' @param gt The ground truth.
#' @param norm Whether or not the confusion matrix is normalized
#' @param graph Whether or not a graphic is displayed.
#' @param ... Other parameters, ignored. \code{\link{performance}} forwards to its evaluation
#' function the same \code{...} it forwards to the learning method, so a call such as
#' \code{performance (BAGGING, ..., type = "confusion", learningmethod = CART)} would otherwise
#' fail on the unused \code{learningmethod}.
#' @return The confusion matrix.
#' @export
#' @seealso \code{\link{evaluation}}, \code{\link{performance}}, \code{\link{splitdata}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = splitdata (iris, 5)
#' model = NB (d$train.x, d$train.y)
#' pred = predict (model, d$test.x)
#' confusion (pred, d$test.y)
confusion <-
  function (predictions, gt, norm = TRUE, graph = TRUE, ...)
  {
    # align.labels() keeps the matrix square even when the model never predicts one of the
    # classes -- otherwise table() silently drops the missing column and the row/column
    # indices below stop referring to the same classes.
    a = align.labels (predictions, gt)
    predictions = a$predictions
    gt = a$gt
    conf = table (gt, predictions, dnn = c ("True labels", "Predicted labels"))
    color = NULL
    maxval = 1
    if (norm)
    {
      rs = rowSums (conf)
      rs [rs == 0] = 1 # a class that never occurs in the ground truth would give 0/0 = NaN
      conf = sweep (conf, 1, rs, "/")
      color = (conf * 100) + 1
    } else {
      color = round (100 * conf / max (conf)) + 1
      maxval = max (conf)
    }
    if (graph)
    {
      graphics::layout (matrix (1:2, ncol = 2), width = c (2, 1), height = c (1, 1))
      on.exit (graphics::layout (1))
      palette = grDevices::colorRampPalette (c ("#FAFAFF", "blue")) (101)
      graphics::plot (c (0, ncol (conf) + 1), c (0, nrow (conf) + 1), col = 0,
                      xlim = c (1, ncol (conf) + 1), ylim = c (1, nrow (conf) + 1),
                      xaxs = "i", yaxs = "i", xlab = "", ylab = "",
                      asp = 1, axes = FALSE,
                      main = "Confusion matrix")
      for (rrow in 1:nrow (conf))
        for (col in 1:ncol (conf))
        {
          row = 1 + nrow (conf) - rrow
          graphics::polygon (x = c (col, col, col + 1, col + 1), y = c (row, row + 1, row + 1, row),
                             col = palette [color [rrow, col]], border = FALSE)
          graphics::text (col + .5, row + .5, round (conf [rrow, col], 2))
        }
      cex = min (5, .15 / log10 (1 + max (graphics::strwidth (c (levels (gt), levels (predictions))))))
      graphics::polygon (x = c (1, 1, ncol (conf) + 1, ncol (conf) + 1), y = c (1, nrow (conf) + 1, nrow (conf) + 1, 1))
      # The horizontal axis carries the columns of 'conf' (the predicted labels) and the
      # vertical one its rows (the true labels), bottom-up hence the rev(). Reading the names
      # off the table itself, rather than off levels(gt) / levels(predictions), keeps the two
      # axes from being swapped.
      graphics::axis (side = 1, at = seq (1.5, by = 1, length.out = ncol (conf)), lwd = 0, lwd.ticks = 1,
                      labels = colnames (conf), pos = 1, cex.axis = cex)
      graphics::axis (side = 2, at = seq (1.5, by = 1, length.out = nrow (conf)), lwd = 0, lwd.ticks = 1,
                      labels = rev (rownames (conf)), pos = 1, cex.axis = cex)
      graphics::title (xlab = "Predicted labels")
      graphics::mtext ("True labels", side = 2, line = 3)

      raster = grDevices::as.raster (matrix (rev (palette), ncol = 1))
      graphics::plot (c (0, 3), c (0, 1), type = 'n', axes = F, xlab = '', ylab = '')
      graphics::rasterImage (raster, 0, 0, 1, 1)
      labels = NULL
      if (maxval == 1)
        labels = seq (0, 1, l = 5)
      else
        labels = round (seq (0, maxval, l = 5))
      graphics::axis (side = 4, at = seq (0, 1, l = 5), lwd = 0, lwd.ticks = 1,
                      labels = labels, pos = 1, cex.axis = 1, las=2)
    }
    return (conf)
  }

#' @keywords internal
# Whether a numeric matrix holds nothing but the integer level codes of 'gt' -- the residue of
# cbind() on factors. Genuine scores are probabilities in [0, 1], so a matrix whose values are
# whole numbers drawn from 1:nlevels(gt), with at least one above 1, cannot be one.
islevelcodes <-
  function (m, gt)
  {
    k = nlevels (factor (gt))
    v = as.vector (m)
    v = v [is.finite (v)]
    return ((length (v) > 0) && all (v == round (v)) && all (v >= 1) && all (v <= k) &&
            any (v > 1))
  }

#' @keywords internal
# Turns whatever the caller passed as 'predictions' into a list of score vectors, one per
# method, oriented so that a high score means "positive class" -- what ROCR expects.
#
# Four shapes are recognised:
#   * a factor or a character vector : hard labels, one method
#   * a numeric vector               : scores of the positive class, one method
#   * a matrix or data.frame whose column names are the class labels -- what
#     predict (model, x, fuzzy = TRUE) returns : the column of 'positive' is used
#   * any other matrix or data.frame : one column per method, each of the two kinds above
#
# 'type' forces the reading: "hard" reduces probabilities to the predicted class first,
# "fuzzy" refuses hard labels, "auto" (the default) goes by shape.
curve.scores <-
  function (predictions, gt, positive, type = c ("auto", "fuzzy", "hard"))
  {
    type = match.arg (type [1], c ("auto", "fuzzy", "hard"))
    ishard = function (v) is.factor (v) || is.character (v) || is.logical (v)
    score = function (v, name)
    {
      if (ishard (v))
      {
        if (type == "fuzzy")
          stop ("roc.curves/cost.curves: type = \"fuzzy\" needs scores (the estimated ",
                "probability of the positive class, e.g. predict (model, x, fuzzy = TRUE)), ",
                "but ", name, " holds hard class labels. Use type = \"hard\" to get the ",
                "(deliberately coarse) curve such labels can produce.")
        return (as.numeric (as.character (v) == positive))
      }
      return (as.numeric (v))
    }
    if (is.null (dim (predictions)))
      return (list (score (predictions, "'predictions'")))
    # cbind() of two factors -- the documented idiom for comparing two methods' hard labels --
    # drops the factor class and leaves the integer level codes. Read as scores, those codes
    # rank the observations by level number rather than by the positive class, so the curve
    # came out upside down whenever 'positive' was not the last level, and type = "hard"
    # compared "1"/"2" against the class labels and scored everything 0. Map them back.
    if (is.numeric (predictions) && islevelcodes (predictions, gt))
      predictions = matrix (levels (factor (gt)) [predictions], nrow = nrow (predictions),
                            dimnames = dimnames (predictions))
    cn = colnames (predictions)
    if ((!is.null (cn)) && (all (levels (factor (gt)) %in% cn)))
    {
      # A matrix of class probabilities, i.e. one column per class, not per method.
      if (type == "hard")
        return (list (as.numeric (cn [apply (predictions, 1, which.max)] == positive)))
      return (list (as.numeric (predictions [, positive])))
    }
    if (type == "hard")
      return (lapply (seq_len (ncol (predictions)),
                      function (j) score (as.character (predictions [, j]),
                                          paste ("column", j, "of 'predictions'"))))
    return (lapply (seq_len (ncol (predictions)),
                    function (j) score (predictions [, j],
                                        paste ("column", j, "of 'predictions'"))))
  }

#' @keywords internal
# Shared plotting loop for roc.curves() and cost.curves().
curve.plot <-
  function (predictions, gt, methods.names, positive, type, measure, ...)
  {
    gt = factor (gt)
    if (nlevels (gt) != 2)
      stop ("roc.curves/cost.curves only work on two-class problems; 'gt' has ",
            nlevels (gt), " classes.")
    if (!(positive %in% levels (gt)))
      stop ("roc.curves/cost.curves: 'positive' (\"", positive, "\") is not one of the class ",
            "labels (", paste (levels (gt), collapse = ", "), ").")
    scores = curve.scores (predictions, gt, positive, type)
    labels = as.numeric (gt == positive)
    for (i in seq_along (scores))
    {
      pred = ROCR::prediction (scores [[i]], labels)
      perf = if (is.na (measure [2])) ROCR::performance (pred, measure [1])
             else ROCR::performance (pred, measure [1], measure [2])
      if (i == 1)
        ROCR::plot (perf, ...)
      else
        ROCR::plot (perf, add = TRUE, lty = i, col = i)
    }
    if ((length (scores) > 1) && (!is.null (methods.names)))
      graphics::legend (if (measure [1] == "tpr") "bottomright" else "topleft",
                        methods.names, lty = seq_along (scores), col = seq_along (scores),
                        bty = "n")
    invisible (NULL)
  }

#' Plot Cost Curves
#'
#' This function plots Cost Curves of several classification predictions.
#' @name cost.curves
#' @inheritParams roc.curves
#' @return Nothing; the curves are drawn on the current graphics device.
#' @export
#' @seealso \code{\link{roc.curves}}, \code{\link{performance}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' model.nb = NB (d [, -5], d [, 5])
#' model.lda = LDA (d [, -5], d [, 5])
#' # From the estimated probabilities (the meaningful version)
#' cost.curves (predict (model.nb, d [, -5], fuzzy = TRUE), d [, 5])
#' # From hard labels, for comparison
#' cost.curves (cbind (predict (model.nb, d [, -5]), predict (model.lda, d [, -5])),
#'              d [, 5], c ("NB", "LDA"), type = "hard")
cost.curves <-
  function (predictions, gt, methods.names = NULL, positive = levels (factor (gt)) [1],
            type = c ("auto", "fuzzy", "hard"), ...)
  {
    curve.plot (predictions, gt, methods.names, positive, type,
                measure = c ("ecost", NA), xlab = "", ylab = "Error", ...)
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
# Puts a vector of predictions and a vector of ground truth values on a *common* set of
# factor levels, so that table(), diag() and label-based indexing behave consistently even
# when one of the two never takes some of the values of the other.
#
# This matters more than it looks: table() only keeps the values it actually observes, so
# table (predictions, gt) is not square as soon as a model never predicts one of the classes
# (routine on imbalanced data, and systematic on small cross-validation folds). Every
# position-based computation on such a table -- diag() above all -- then reads cells that
# have nothing to do with the intended ones, and returns a plausible but wrong number
# instead of raising an error.
align.labels <-
  function (predictions, gt)
  {
    if (is.factor (predictions) || is.factor (gt) || is.character (predictions) || is.character (gt))
      lev = union (levels (factor (gt)), levels (factor (predictions)))
    else
      lev = sort (union (unique (gt), unique (predictions)))
    return (list (predictions = factor (predictions, levels = lev),
                  gt = factor (gt, levels = lev),
                  levels = lev))
  }

#' @keywords internal
# Shared preparation for LDA() and QDA(): checks that every predictor is numeric, and returns
# the logical vector of the columns with a non-zero variance -- the ones the model can
# actually be fitted on.
#
# Indexing with that vector must always use drop = FALSE: without it, a single retained
# predictor collapses to a bare vector, which MASS::lda()/qda() and their predict methods
# reject with "'x' is not a matrix" -- and a discriminant analysis on one variable is a
# textbook exercise.
discriminant.variables <-
  function (train, fname)
  {
    d = as.data.frame (train)
    numeric = sapply (d, is.numeric)
    if (!all (numeric))
      stop (fname, ": every predictor must be numeric -- discriminant analysis cannot use ",
            "qualitative variables. Non-numeric column(s): ",
            paste (colnames (d) [!numeric], collapse = ", "))
    variables = apply (train, 2, stats::var) > 0
    if (!any (variables))
      stop (fname, ": every predictor is constant (zero variance), there is nothing to ",
            "discriminate on.")
    return (variables)
  }

#' @keywords internal
# Per-class counts behind precision, recall and everything derived from them. Rows of the
# table are the ground truth, columns the predictions, on the union of the two sets of labels
# (see align.labels): 'tp' is its diagonal, 'predicted' the column totals (TP + FP) and
# 'actual' the row totals (TP + FN).
prf.counts <-
  function (predictions, gt)
  {
    a = align.labels (predictions, gt)
    t = table (a$gt, a$predictions)
    return (list (tp = diag (t), predicted = colSums (t), actual = rowSums (t),
                  levels = a$levels, classes = levels (factor (gt)), n = sum (t)))
  }

#' @keywords internal
# x / y, with the 0 / 0 that a class nobody predicts (or nobody belongs to) produces read as
# 0 rather than NaN -- the usual convention, and the one that keeps an average finite.
prf.ratio <-
  function (x, y) ifelse (y == 0, 0, x / y)

#' @keywords internal
# Every measure in this family is a function of a precision and a recall; this is the only
# place each formula is written.
prf.value <-
  function (precision, recall, what, beta = 1)
  {
    if (what == "precision")
      return (precision)
    if (what == "recall")
      return (recall)
    if (what == "fmeasure")
      return (prf.ratio ((1 + beta * beta) * precision * recall,
                         beta * beta * precision + recall))
    if (what == "goodness")
      return ((beta * precision + recall) / (beta + 1))
    if (what == "jaccard")
      return (prf.ratio (precision * recall, precision + recall - precision * recall))
    if (what == "fowlkesmallows")
      return (sqrt (precision * recall))
    stop ("prf.value: unknown measure '", what, "'.")
  }

#' @keywords internal
# Shared body of evaluation.precision(), evaluation.recall(), evaluation.fmeasure(),
# evaluation.goodness(), evaluation.jaccard() and evaluation.fowlkesmallows().
#
# These measures are defined for one class against the rest. On a two-class problem that class
# is 'positive' and there is nothing to average; beyond two classes the per-class values have
# to be combined, which is what 'average' selects. Averages run over the classes of the ground
# truth: a label that only ever appears among the predictions is a mistake, not a class of the
# problem, and it is counted as such (it lowers the recall of the classes it was predicted
# instead of) without contributing a term of its own.
prf <-
  function (predictions, gt, what, average = NULL, positive = NULL, beta = 1)
  {
    gt = factor (gt)
    counts = prf.counts (predictions, gt)
    classes = counts$classes
    if (is.null (average))
      average = if (length (classes) == 2) "binary" else "macro"
    average = match.arg (average [1], c ("binary", "macro", "micro", "weighted", "none"))
    if (average == "binary")
    {
      if (length (classes) != 2)
        stop ("evaluation.", what, ": average = \"binary\" needs exactly two classes, but the ",
              "ground truth has ", length (classes), " (",
              paste (classes, collapse = ", "), "). Use average = \"macro\", \"micro\", ",
              "\"weighted\" or \"none\" -- \"macro\" is the default beyond two classes.")
      if (is.null (positive))
        positive = classes [1]
      if (!(positive %in% classes))
        stop ("evaluation.", what, ": positive = \"", positive, "\" is not one of the class ",
              "labels (", paste (classes, collapse = ", "), ").")
      p = prf.ratio (counts$tp [positive], counts$predicted [positive])
      r = prf.ratio (counts$tp [positive], counts$actual [positive])
      return (unname (prf.value (p, r, what, beta)))
    }
    if (average == "micro")
    {
      # Pooling the counts of every class before dividing. With single-label predictions each
      # observation contributes exactly one predicted and one actual label, so the two totals
      # are both the number of observations: micro-averaged precision, recall, F-measure and
      # accuracy all coincide. Worth showing students rather than hiding.
      p = prf.ratio (sum (counts$tp), counts$n)
      return (unname (prf.value (p, p, what, beta)))
    }
    p = prf.ratio (counts$tp [classes], counts$predicted [classes])
    r = prf.ratio (counts$tp [classes], counts$actual [classes])
    value = prf.value (p, r, what, beta)
    names (value) = classes
    if (average == "none")
      return (value)
    if (average == "macro")
      return (mean (value))
    return (stats::weighted.mean (value, counts$actual [classes]))
  }

#' Shared documentation for the 'average' and 'positive' parameters
#'
#' This function is never called: it holds the canonical documentation of the \code{average}
#' and \code{positive} parameters, shared (via \code{@@inheritParams}) by the six evaluation
#' measures built on a precision and a recall.
#' @name average.doc
#' @param average How the per-class values are combined. These measures are defined for one
#' class against all the others, so a single number requires either picking that class or
#' averaging.
#' \describe{
#'   \item{\code{"binary"}}{the value for the class named by \code{positive}. Two-class
#'     problems only, and the default there.}
#'   \item{\code{"macro"}}{the plain mean of the per-class values. The default beyond two
#'     classes. Gives every class the same weight, whatever its size, so a rare class the model
#'     never gets right weighs as much as the majority one.}
#'   \item{\code{"weighted"}}{the mean of the per-class values, weighted by the number of
#'     observations of each class.}
#'   \item{\code{"micro"}}{pools the counts of every class before dividing. With single-label
#'     predictions each observation contributes one predicted and one actual label, so
#'     micro-averaged precision, recall and F-measure all equal the accuracy.}
#'   \item{\code{"none"}}{the vector of per-class values, named after the classes. Cannot be
#'     used through \code{\link{evaluation}} or \code{\link{performance}}, which expect one
#'     number per criterion.}
#' }
#' Averages run over the classes of the ground truth. A class that no observation is predicted
#' to belong to has an undefined precision; it is read as 0, the usual convention.
#' @param positive The label of the positive class, used by \code{average = "binary"} only.
#' Defaults to \code{levels (gt) [1]}, i.e. the \emph{first} level of the ground truth factor
#' -- which, for the usual alphabetical level ordering, is often the \emph{negative} class
#' (\code{"N"} before \code{"Y"}, \code{"No"} before \code{"Yes"}, ...). Set this argument
#' explicitly whenever the positive class is not the first level.
#' @keywords internal
average.doc <-
  function (average, positive)
  {
    NULL
  }

#' @keywords internal
eval.accuracy <-
  function (predictions, gt, ...)
  {
    a = align.labels (predictions, gt)
    return (sum (a$predictions == a$gt, na.rm = TRUE) / length (a$gt))
  }

#' @keywords internal
eval.fmeasure <-
  function (predictions, gt, beta = 1, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "fmeasure", average, positive, beta)

#' @keywords internal
eval.fowlkesmallows <-
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "fowlkesmallows", average, positive)

#' @keywords internal
eval.goodness <-
  function (predictions, gt, beta = 1, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "goodness", average, positive, beta)

#' @keywords internal
eval.jaccard <-
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "jaccard", average, positive)

#' @keywords internal
# Cohen's kappa: the agreement between the predictions and the ground truth, corrected for the
# agreement two independent labellings would reach by chance.
#
# The two vectors go through align.labels() first, so that a model which never predicts some
# class is not compared on shifted codes. Class labels are nominal, so this is the *unweighted*
# kappa: every mistake counts the same, whichever two classes were confused.
eval.kappa <-
  function (predictions, gt, ...)
  {
    a = align.labels (predictions, gt)
    t = table (a$gt, a$predictions)
    n = sum (t)
    if (n == 0)
      return (NA_real_)
    po = sum (diag (t)) / n
    pe = sum (rowSums (t) * colSums (t)) / (n * n)
    # A single class, always predicted: chance agreement is already perfect, and kappa is
    # 0 / 0 rather than 1.
    if (pe >= 1)
      return (NA_real_)
    return ((po - pe) / (1 - pe))
  }

#' @keywords internal
eval.precision <-
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "precision", average, positive)

#' @keywords internal
eval.recall <-
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "recall", average, positive)

#' @keywords internal
# Lookup table for the eval.* dispatch used by evaluation(). Built lazily (as a function
# rather than a static list) so it does not depend on file collation order: eval.adjr2,
# eval.msep and eval.r2 are defined in regression.R, the others in this file.
eval.functions <-
  function ()
  {
    list (accuracy = eval.accuracy,
          adjr2 = eval.adjr2,
          fmeasure = eval.fmeasure,
          fowlkesmallows = eval.fowlkesmallows,
          goodness = eval.goodness,
          jaccard = eval.jaccard,
          kappa = eval.kappa,
          msep = eval.msep,
          precision = eval.precision,
          r2 = eval.r2,
          recall = eval.recall)
  }

#' Evaluation of classification or regression predictions
#'
#' Evaluation predictions of a classification or a regression model.
#' @name evaluation
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth of the dataset (\code{factor} or \code{vector}).
#' @param eval The evaluation method.
#' @param ... Other parameters.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{confusion}}, \code{\link{evaluation.accuracy}}, \code{\link{evaluation.fmeasure}}, \code{\link{evaluation.fowlkesmallows}}, \code{\link{evaluation.goodness}}, \code{\link{evaluation.jaccard}}, \code{\link{evaluation.kappa}},
#' \code{\link{evaluation.precision}}, \code{\link{evaluation.recall}},
#' \code{\link{evaluation.msep}}, \code{\link{evaluation.r2}}, \code{\link{performance}}
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
  function (predictions, gt, eval = ifelse (is.factor (gt), "accuracy", "r2"), ...)
  {
    # Each criterion computes what it needs from 'predictions' and 'gt' -- one contingency
    # table each. A macro-averaged measure is an average over the classes, not a measure
    # derived from aggregated counts, so it cannot be built from a shared precision and recall.
    funs = eval.functions ()
    res = NULL
    for (e in eval)
    {
      fun = funs [[e]]
      if (is.null (fun))
        stop ("evaluation: unknown evaluation criterion '", e, "'. Available criteria: ",
              paste (names (funs), collapse = ", "))
      tmp = fun (predictions = predictions, gt = gt, ...)
      if (length (tmp) != 1)
        stop ("evaluation: criterion '", e, "' returned ", length (tmp), " values instead of ",
              "one. If this is average = \"none\", call evaluation.", e, "() directly: ",
              "evaluation() returns one number per criterion.")
      res = c (res, tmp)
    }
    names (res) = eval
    return (res)
  }

#' Accuracy of classification predictions
#'
#' Evaluation predictions of a classification model according to accuracy.
#' @name evaluation.accuracy
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
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
  function (predictions, gt, ...)
  {
    return (eval.accuracy (predictions, gt))
  }

#' F-measure
#'
#' Evaluation predictions of a classification model according to the F-measure index.
#' @name evaluation.fmeasure
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @param beta The weight given to precision.
#' @inheritParams average.doc
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
  function (predictions, gt, beta = 1, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "fmeasure", average, positive, beta)

#' Fowlkes–Mallows index
#'
#' Evaluation predictions of a classification model according to the Fowlkes–Mallows index.
#' @name evaluation.fowlkesmallows
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @inheritParams average.doc
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
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "fowlkesmallows", average, positive)

#' Goodness
#'
#' Evaluation predictions of a classification model according to Goodness index.
#' @name evaluation.goodness
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @param beta The weight given to precision.
#' @inheritParams average.doc
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
  function (predictions, gt, beta = 1, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "goodness", average, positive, beta)

#' Jaccard index
#'
#' Evaluation predictions of a classification model according to Jaccard index.
#' @name evaluation.jaccard
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @inheritParams average.doc
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
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "jaccard", average, positive)

#' Kappa evaluation of classification predictions
#'
#' Evaluation predictions of a classification model according to Cohen's kappa: the proportion
#' of correct predictions, corrected for the proportion two independent labellings would get
#' right by chance. Class labels being nominal, the kappa is the unweighted one -- every
#' mistake counts the same.
#' @name evaluation.kappa
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
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
  function (predictions, gt, ...)
  {
    return (eval.kappa (predictions, gt))
  }

#' Precision of classification predictions
#'
#' Evaluation predictions of a classification model according to precision.
#' @name evaluation.precision
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @inheritParams average.doc
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
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "precision", average, positive)

#' Recall of classification predictions
#'
#' Evaluation predictions of a classification model according to recall.
#' @name evaluation.recall
#' @param predictions The predictions of a classification model (\code{factor} or \code{vector}).
#' @param gt The ground truth (\code{factor} or \code{vector}).
#' @inheritParams average.doc
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
  function (predictions, gt, average = NULL, positive = NULL, ...)
    prf (predictions, gt, "recall", average, positive)

#' Classification using Gradient Boosting
#'
#' This function builds a classification model using Gradient Boosting
#' @name GRADIENTBOOSTING
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param ntree The number of trees in the forest.
#' @param learningrate The learning rate (between 0 and 1).
#' @param seed A specified seed for random number generation (row/column subsampling, if used
#' via \code{...}; \code{xgboost}'s default parameters are otherwise deterministic, but the seed
#' is provided for consistency with the rest of the package's API and to cover subsampling
#' parameters passed through \code{...}).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{GRADIENTBOOSTING}
#' does not yet implement hyperparameter tuning.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{GRADIENTBOOSTING} does not produce
#' a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[xgboost]{xgboost}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' GRADIENTBOOSTING (iris [, -5], iris [, 5])
#' }
GRADIENTBOOSTING <-
  function (train, labels,
            ntree = 500,
            learningrate = 0.3,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      # Since xgboost >= 2.1, the data/label/eta/verbose parameter names are deprecated in favour
      # of x/y/learning_rate/verbosity, and the classification objective (binary:logistic or
      # multi:softprob) is inferred automatically from 'y' being a factor: passing a 0/1-recoded
      # numeric label (the previous approach) now errors out instead of being auto-detected.
      setseed (seed)
      model = xgboost::xgboost (x = as.matrix (train), y = labels, nrounds = ntree, learning_rate = learningrate, verbosity = 0)
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
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{KNN} does not
#' support reusing pre-tuned parameters (it stores the training set and re-tunes \code{k} on
#' every call when \code{k} is a vector).
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{KNN} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[class]{knn}}
#' @examples
#' require (datasets)
#' data (iris)
#' KNN (iris [, -5], iris [, 5])
KNN <-
  function (train, labels, k = 1:10, nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      setseed (seed)
      kk = k [1]
      if (is.vector (k) && (length (k) > 1))
      {
        tunecontrol = tune.scheme (nfolds)
        # '$best.parameters$k', not '$best.model$k': e1071 documents the former, and the
        # latter is empty whenever tune.knn() could not refit a best model. NULL then went
        # into the result, where list() silently *drops* it -- so the model had no 'k' at all
        # and predict() failed much later on "k = 0 must be at least 1".
        tuned = e1071::tune.knn (train, labels, k = k, tunecontrol = tunecontrol)
        kk = tuned$best.parameters$k
        if (is.null (kk) || is.na (kk) || (length (kk) != 1) || (kk < 1))
        {
          kk = k [1]
          warning ("KNN: the cross-validated choice of k did not return a usable value; ",
                   "falling back on k = ", kk, ".")
        }
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
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{LDA} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{LDA} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[MASS]{lda}}
#' @examples
#' require (datasets)
#' data (iris)
#' LDA (iris [, -5], iris [, 5])
LDA <-
  function (train, labels, tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      variables = discriminant.variables (train, "LDA")
      model = MASS::lda (x = train [, variables, drop = FALSE], grouping = labels)
      res = list (model = model, method = "LDA", variables = variables)
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
#' @param reg The penalty applied to the coefficients, as in \code{\link{LINREG}}:
#' \code{"none"} (the default) fits the plain multinomial logistic regression of
#' \code{\link[nnet]{multinom}}, while \code{"ridge"} (L2), \code{"lasso"} (L1) and
#' \code{"elastic"} (a mix of the two, weighted by \code{alpha}) fit a penalized one with
#' \code{\link[glmnet]{glmnet}}. Penalizing is what makes logistic regression usable when the
#' predictors are numerous or strongly correlated, where the unpenalized fit either fails to
#' converge or separates the classes perfectly with unbounded coefficients.
#' @param lambda The grid of penalty strengths searched by cross-validation; the retained value
#' is the one minimising the cross-validated deviance. \code{NULL} (the default) lets
#' \code{\link[glmnet]{glmnet}} derive the grid from the data, which is the recommended
#' choice: a fixed grid reaching very small penalties makes the fit fail to converge on
#' separable data.
#' @param alpha The elastic net mixing parameter, between 0 (ridge) and 1 (lasso). Used by
#' \code{reg = "elastic"} only.
#' @param nfolds The number of folds of the cross-validation used to choose \code{lambda}.
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{LR} does not
#' support reusing pre-tuned parameters.
#' @param graph Whether the cross-validation curve used to choose \code{lambda} is plotted.
#' Ignored by \code{reg = "none"}, which has nothing to choose.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[nnet]{multinom}}, \code{\link[glmnet]{glmnet}}, \code{\link{LINREG}}
#' @examples
#' require (datasets)
#' data (iris)
#' LR (iris [, -5], iris [, 5])
#' \donttest{
#' # Penalized variants: same three penalties as LINREG()
#' d = splitdata (iris, 5, seed = 0)
#' model = LR (d$train.x, d$train.y, reg = "lasso")
#' evaluation (predict (model, d$test.x), d$test.y)
#' }
LR <-
  function (train, labels, reg = c ("none", "ridge", "lasso", "elastic"),
            lambda = NULL, alpha = .5, nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      # See CDA: the previous guard rejected character labels outright (nlevels() == 0) and
      # returned NULL silently.
      labels = check.classes (labels, "LR")
      if (is.vector (train))
      {
        train = matrix (train, ncol = 1)
        colnames (train) = "X"
      }
      reg = match.arg (reg [1], c ("none", "ridge", "lasso", "elastic"))
      if (reg == "none")
      {
        data = cbind.data.frame (train, Class = labels)
        model = nnet::multinom (formula = Class~., data, trace = FALSE)
        res = list (model = model, method = "LR")
      }
      else
      {
        # Penalized logistic regression, the classification counterpart of
        # LINREG (reg = "ridge" / "lasso" / "elastic"): same three penalties, same meaning of
        # 'alpha', and lambda chosen by cross-validation on the same grid.
        palpha = switch (reg, ridge = 0, lasso = 1, alpha)
        # Left to glmnet by default: it derives the grid from the data, starting at the
        # smallest penalty that shrinks every coefficient to zero. A fixed grid reaching very
        # small penalties (the one LINREG() uses) makes the logistic fit fail to converge on
        # separable data, which is exactly the situation penalization is meant to handle, and
        # fills the console with warnings.
        if (!is.null (lambda))
          lambda = sort (lambda, decreasing = TRUE)
        family = if (nlevels (labels) == 2) "binomial" else "multinomial"
        # glmnet warns, once per offending value of the grid, that the fit did not converge
        # for the smallest penalties. On separable data that is expected rather than wrong --
        # unpenalized logistic regression has no finite solution there, which is the very
        # reason for penalizing -- and cv.glmnet() picks its lambda among the values that did
        # converge. The raw C++ warnings are replaced by one sentence saying so.
        notconverged = 0
        cv = withCallingHandlers (
          glmnet::cv.glmnet (as.matrix (train), labels, family = family, alpha = palpha,
                             lambda = lambda, standardize = TRUE, nfolds = nfolds),
          warning = function (w)
          {
            if (grepl ("Convergence for", conditionMessage (w), fixed = TRUE))
            {
              notconverged <<- notconverged + 1
              invokeRestart ("muffleWarning")
            }
          })
        if (notconverged > 0)
          message ("LR: the fit did not converge for the ", notconverged, " smallest penalties ",
                   "of the grid, which happens when the classes are separable -- an ",
                   "unpenalized logistic regression has no finite solution then. The penalty ",
                   "is chosen among the values that did converge; no action is needed.")
        if (graph)
          graphics::plot (cv)
        model = glmnet::glmnet (as.matrix (train), labels, family = family, alpha = palpha,
                                lambda = cv$lambda.min, standardize = TRUE)
        res = list (model = model, lev = levels (labels), lambda = cv$lambda.min,
                    alpha = palpha, method = "penalizedlr")
      }
      class (res) = "model"
    }
    return (res)
  }

#' Classification using Multilayer Perceptron
#'
#' This function builds a classification model using Multilayer Perceptron.
#' @name MLP
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param hidden The size of the hidden layer (if a vector, cross-over validation is used to chose the best size).
#' @param decay The decay (between 0 and 1) of the backpropagation algorithm (if a vector, cross-over validation is used to chose the best size).
#' @param methodparameters Object containing the parameters. If given, it replaces \code{size} and \code{decay}.
#' @inheritParams tune.doc
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{MLP} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[nnet]{nnet}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' MLP (iris [, -5], iris [, 5], hidden = 4, decay = .1)
#' }
MLP <-
  function (train,
            labels,
            # if/else, not ifelse(): ifelse() is vectorised and would return only the
            # *first* element of the branch it selects, collapsing this grid to one value.
            hidden = if (is.vector (train)) 2:(1 + nlevels (labels)) else 2:(ncol (train) + nlevels (labels)),
            decay = 10^(-3:-1),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    model = NULL
    if (is.vector (train))
      train = data.frame (X = train)
    d = cbind.data.frame (Class = labels, train)
    if (!is.null (methodparameters))
    {
      # Each value is taken only if the 'params' object actually carries it: an *empty* one
      # -- what a method with nothing to tune returns, and what performance() then hands back
      # to it -- must leave the defaults alone rather than overwrite them with NULL.
      if (!is.null (methodparameters$hidden))
        hidden = methodparameters$hidden
      if (!is.null (methodparameters$decay))
        decay = methodparameters$decay
    }
    if (length (hidden) > 1 | length (decay) > 1)
    {
      model = e1071::tune.nnet (Class~., data = d, size = hidden, decay = decay,
                                tunecontrol = tune.scheme (nfolds), ...)$best.model
    }
    else
      model = nnet::nnet (Class~., data = d, size = hidden, decay = decay, trace = FALSE, ...)
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
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{NB} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{NB} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{naiveBayes}}
#' @examples
#' require (datasets)
#' data (iris)
#' NB (iris [, -5], iris [, 5])
NB <-
  function (train, labels, tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
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

#' @keywords internal
panel.compare <-
  function (x, y, ...)
  {
    graphics::points (x, y, col = "red")
    graphics::abline (a = 0, b = 1, col = "blue")
  }

#' Performance estimation
#'
#' Estimate the performance of classification or regression methods using bootstrap or crossvalidation (accuracy, ROC curves, confusion matrices, ...)
#' @name performance
#' @param methods The classification or regression methods to be evaluated.
#' @param train.x The dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' @param train.y The target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param test.x The test dataset (description/predictors), a \code{matrix} or \code{data.frame}.
#' Giving \code{test.x} and \code{test.y} is the simplest use of this function: each method is
#' fitted on \code{(train.x, train.y)}, used to predict \code{test.x}, and its predictions are
#' compared with \code{test.y} -- no resampling at all. \code{protocol} then defaults to
#' \code{"holdout"}, and any other protocol raises an error, since it would ignore the test set.
#' @param test.y The (test) target (class labels or numeric values), a \code{factor} or \code{vector}.
#' @param train.size The size of the training set, for \code{protocol = "holdout"} without an
#' explicit test set: either a number of observations, or a proportion between 0 and 1.
#' @param type The type of evaluation (confusion matrix, ROC curve, ...)
#' @param protocol How the performance is estimated.
#' \describe{
#'   \item{\code{"bootstrap"}}{(default) \code{nruns} bootstrap samples of the training set,
#'     each model evaluated on the observations left out of its sample.}
#'   \item{\code{"crossvalidation"}}{\code{nruns} repetitions of a \code{nfolds}-fold
#'     cross-validation of the training set.}
#'   \item{\code{"loocv"}}{leave-one-out cross-validation.}
#'   \item{\code{"holdout"}}{a single train/test split. Uses \code{test.x}/\code{test.y} when
#'     they are given -- see \code{test.x} -- and otherwise draws a training set of
#'     \code{train.size} observations and evaluates on the rest.}
#'   \item{\code{"train"}}{evaluates each model on the very data it was fitted on. Optimistic
#'     by construction; useful to show students exactly that.}
#' }
#' @param eval The evaluation functions.
#' @param nruns The number of bootstrap runs.
#' @param nfolds The number of folds (crossvalidation estimation).
#' @param new A logical value indicating whether a new plot should be created or not (cost curves or ROC curves).
#' @param lty The line type (and color) specified as an integer (cost curves or ROC curves).
#' @param methodparameters Method parameters (if null tuning is done by cross-validation).
#' @param names Method names.
#' @param fuzzy Used by \code{type = "roc"} and \code{type = "cost"} only. \code{FALSE} by
#' default: the curves are built from the hard class labels, which reduces them to three
#' points. Pass \code{fuzzy = TRUE} to build them from the estimated probabilities of the
#' positive class, which is what a ROC curve is meant to show. See \code{\link{roc.curves}}.
#' @param positive The label of the positive class. Used by \code{type = "roc"} and
#' \code{type = "cost"} to orient the curves, and passed on to \code{\link{evaluation}} for the
#' criteria that are defined on one class -- precision, recall, the F-measure and the other
#' measures taking an \code{average} argument -- so that a two-class problem can be scored on
#' either of its classes. Defaults to the first level of the target, which is worth setting
#' explicitly whenever the class of interest is not the first one.
#' @param stratify Whether the splits should preserve the proportions of the classes
#' (\code{TRUE}, the default), for \code{protocol = "crossvalidation"} and for
#' \code{protocol = "holdout"} when it draws its own split. Ignored for a numeric target, and
#' for the bootstrap and leave-one-out protocols, which have no split to stratify.
#' @param seed A specified seed for random number generation (useful for testing different method with the same bootstap samplings).
#' @param ... Other specific parameters for the leaning method.
#' @return The evaluation of the predictions (numeric value).
#' @export
#' @seealso \code{\link{confusion}}, \code{\link{evaluation}}, \code{\link{cost.curves}}, \code{\link{roc.curves}}
#' @examples
#' \dontrun{
#' require ("datasets")
#' data (iris)
#' # The simplest use: a training set, a test set, and the score of the model fitted on the
#' # first and evaluated on the second. Same thing as
#' # evaluation.accuracy (predict (NB (d$train.x, d$train.y), d$test.x), d$test.y).
#' d = splitdata (iris, 5, seed = 0)
#' performance (NB, d$train.x, d$train.y, d$test.x, d$test.y)
#' # Several methods and criteria at once
#' performance (c (NB, LDA, CART), d$train.x, d$train.y, d$test.x, d$test.y,
#'              eval = c ("accuracy", "kappa"))
#' # One method, one evaluation criterion, bootstrap estimation
#' performance (NB, iris [, -5], iris [, 5], seed = 0)
#' # One method, two evaluation criteria, train set estimation
#' performance (NB, iris [, -5], iris [, 5], eval = c ("accuracy", "kappa"),
#'              protocol = "train", seed = 0)
#' # Three methods, ROC curves, LOOCV estimation
#' data (linsep)
#' performance (c (NB, LDA, LR), linsep [, -3], linsep [, 3], type = "roc",
#'              protocol = "loocv", seed = 0)
#' # Same curves, read from the hard predicted labels instead of the class-membership
#' # scores: each method collapses to a single operating point.
#' performance (c (NB, LDA, LR), linsep [, -3], linsep [, 3], type = "roc",
#'              protocol = "loocv", seed = 0, fuzzy = FALSE)
#' # Choosing the positive class explicitly
#' performance (NB, linsep [, -3], linsep [, 3], type = "roc", protocol = "loocv",
#'              seed = 0, positive = levels (linsep [, 3]) [2])
#' # List of methods in a variable, confusion matrix, hodout estimation
#' classif = c (NB, LDA, LR)
#' performance (classif, iris [, -5], iris [, 5], type = "confusion",
#'              protocol = "holdout", seed = 0, names = c ("NB", "LDA", "LR"))
#' # List of strings (method names), scatterplot evaluation, crossvalidation estimation
#' classif = c ("NB", "LDA", "LR")
#' performance (classif, iris [, -5], iris [, 5], type = "scatter",
#'              protocol = "crossvalidation", seed = 0)
#' # Actual vs. predicted
#' data (trees)
#' performance (LINREG, trees [, -3], trees [, 3], type = "avsp")
#' }
performance <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL, train.size = round (0.7 * nrow (train.x)), type = c ("evaluation", "confusion", "roc", "cost", "scatter", "avsp"),
            protocol = c ("bootstrap", "crossvalidation", "loocv", "holdout", "train"),
            eval = ifelse (is.factor (train.y), "accuracy", "r2"),
            nruns = 10, nfolds = 10, new = TRUE, lty = 1,
            seed = NULL, methodparameters = NULL, names = NULL,
            fuzzy = FALSE, positive = NULL, stratify = TRUE, ...)
  {
    if (length (train.y) == 1)
    {
      index = train.y
      train.y = train.x [, index]
      train.x = train.x [, -index]
    }
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
        # The names of the functions as they were written in the call. A list held in a
        # variable -- 'classif = c (NB, LDA); performance (classif, ...)' -- carries none, and
        # the rows of the result used to come back unnamed; they are numbered instead.
        methodNames = as.character (match.call ()$methods)
        if (length (methodNames) > 1)
          methodNames = methodNames [-1]
        if (length (methodNames) != length (methods))
          methodNames = paste ("Method", seq_along (methods))
      }
    }
    if (is.vector (train.x))
      train.x = data.frame (X = train.x)
    setseed (seed)
    if (is.null (methodparameters))
    {
      if (length (methods) == 1)
        methodparameters = methods (train.x, train.y, tune = TRUE, ...)
      else
        # The seed is set again before each method rather than once for the whole list: a
        # method that searches a grid consumes random numbers, so tuning them in sequence made
        # the hyperparameters retained for a method depend on which methods came before it in
        # the vector. The same method with the same seed then scored differently in
        # 'c (KNN, CART, MLP, SVMr)' and in 'c (SVMl, SVMr)'.
        methodparameters = lapply (methods, function (method)
        {
          setseed (seed)
          method (train.x, train.y, tune = TRUE, ...)
        })
    }
    protocols = protocol.functions ()
    # Only "holdout" uses an explicitly supplied test set; every other protocol builds its own
    # splits out of the training set. So: a test set with no protocol named means evaluate on
    # that test set, and a test set with a protocol that cannot use it is an error rather than
    # a silent no-op.
    hastest = (!is.null (test.x)) && (!is.null (test.y))
    if (hastest)
    {
      if (missing (protocol))
        protocol = "holdout"
      else if (match.arg (protocol [1], names (protocols)) != "holdout")
        stop ("performance: 'test.x' and 'test.y' were given, but protocol = \"", protocol [1],
              "\" builds its own splits from the training set and would ignore them. Use ",
              "protocol = \"holdout\" to fit on (train.x, train.y) and evaluate on ",
              "(test.x, test.y), or drop 'test.x'/'test.y' to estimate the performance by ",
              "resampling the training set alone.")
    }
    protocolfun = protocols [[protocol [1]]]
    if (is.null (protocolfun))
      stop ("performance: unknown protocol '", protocol [1], "'. Available protocols: ",
            paste (names (protocols), collapse = ", "))
    # ROC and cost curves need a score, not a class: ask the protocol for probabilities.
    wantfuzzy = (type [1] %in% c ("roc", "cost")) && fuzzy
    tmp = protocolfun (methods = methods, train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y, train.size = train.size,
                       methodparameters = methodparameters, nruns = nruns, nfolds = nfolds, seed = seed, fuzzy = wantfuzzy,
                       stratify = stratify, ...)
    predictions = tmp$predictions
    targets = tmp$targets
    if (wantfuzzy)
    {
      # predictions is a list over splits of lists over methods; stack the splits, then keep
      # the probability of the positive class for each method.
      if (is.null (positive))
        positive = levels (factor (train.y)) [1]
      if (nlevels (factor (train.y)) != 2)
        stop ("performance: type = \"", type [1], "\" needs a two-class problem; the target has ",
              nlevels (factor (train.y)), " classes.")
      targets = factor (unlist (targets), labels = levels (train.y))
      probs = lapply (seq_len (length (methods)),
                      function (m) do.call (rbind, lapply (predictions, function (split) split [[m]])))
      missing = sapply (probs, function (m) is.null (colnames (m)) || (!(positive %in% colnames (m))))
      if (any (missing))
        stop ("performance: the fuzzy predictions of method ", which (missing) [1],
              " do not carry a column named after the positive class (\"", positive, "\"), ",
              "so no score can be extracted. Use fuzzy = FALSE to fall back on the coarse ",
              "curves built from hard class labels.")
      scores = sapply (probs, function (m) as.numeric (m [, positive]))
      if (length (methods) > 1)
        colnames (scores) = methodNames
      if (type [1] == "roc")
        return (invisible (roc.curves (scores, targets, methodNames, positive = positive, type = "fuzzy")))
      return (invisible (cost.curves (scores, targets, methodNames, positive = positive, type = "fuzzy")))
    }
    if (type [1] %in% c ("evaluation", "confusion", "roc", "cost", "avsp"))
    {
      if (length (methods) == 1)
        predictions = unlist (predictions)
      else
        predictions = do.call ("rbind", predictions)
      predictions = unlist (predictions)
      targets = unlist (targets)
      if (is.factor (train.y))
      {
        lab = levels (train.y)
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
        }
        targets = factor (targets, labels = levels (train.y))
      }
      if (length (methods) > 1)
        predictions = as.data.frame (predictions)
    }
    if (type [1] == "evaluation")
    {
      res = NULL
      if (length (methods) == 1)
        res = evaluation (predictions = predictions, gt = targets, eval = eval, ncol = ncol (train.x), positive = positive, ...)
      else
      {
        res = t (as.data.frame (lapply (predictions, function (column) evaluation (predictions = column, gt = targets, eval = eval, ncol = ncol (train.x), positive = positive, ...))))
        rownames (res) = methodNames
      }
      return (res)
    }
    else if (type [1] == "confusion")
    {
      res = NULL
      if (length (methods) == 1)
        res = confusion (predictions = predictions, gt = targets, ...)
      else
      {
        for (prediction in predictions)
          res = c (res, list (confusion (predictions = prediction, gt = targets, ...)))
        names (res) = methodNames
      }
      return (res)
    }
    else if (type [1] == "avsp")
      plotavsp (predictions, targets)
    else if (type [1] == "roc")
      roc.curves (predictions, targets, methodNames,
                  positive = if (is.null (positive)) levels (factor (targets)) [1] else positive,
                  type = "hard")
    else if (type [1] == "cost")
      cost.curves (predictions, targets, methodNames,
                   positive = if (is.null (positive)) levels (factor (targets)) [1] else positive,
                   type = "hard")
    else if (type [1] == "scatter")
    {
      res = array (dim = c (length (targets), length (methods), length (eval)))
      for (i in 1:length (targets))
      {
        gt = targets [[i]]
        pred = predictions [[i]]
        if (is.factor (train.y))
        {
          lab = levels (train.y)
          l1 = length (lab)
          if (length (methods) == 1)
          {
            l2 = length (unique (pred))
            if (l2 > l1)
              lab = c (lab, rep ("Unknown", l2 - l1))
            lab = lab [sort (unique (pred))]
            pred = factor (pred, labels = lab)
          }
          else
          {
            pred = as.data.frame (pred)
            pred = lapply (pred, function (column) {
              l2 = length (unique (column))
              lab2 = lab
              if (l2 > l1)
                lab2 = c (lab, rep ("Unknown", l2 - l1))
              lab2 = lab2 [sort (unique (column))]
              return (factor (column, labels = lab2))
            })
            pred = as.data.frame (pred)
          }
          gt = factor (gt, labels = levels (train.y))
          if (length (methods) == 1)
            res [i, , ] = as.matrix (evaluation (predictions = pred, gt = gt, eval = eval, ncol = ncol (train.x), positive = positive, ...))
          else
            res [i, , ] = as.matrix (t (as.data.frame (lapply (pred, function (column) evaluation (predictions = column, gt = gt, eval = eval, ncol = ncol (train.x), positive = positive, ...)))))
        }
      }
      for (i in 1:dim (res) [3])
      {
        if (dim (res) [2] == 1)
          graphics::plot (sort (res [, 1, i]), col = "red", xlab = "", ylab = dimnames (res) [[3]][i], main = dimnames (res) [[2]][1])
        else if (dim (res) [2] == 2)
        {
          dd = matrix (res [, 1:2, i], nrow = dim (res) [1])
          colnames (dd) = methodNames
          lim = c (min (dd), max (dd))
          graphics::plot (dd, asp = 1, xlim = lim, ylim = lim, col = "red", main = dimnames (res) [[3]][i])
          graphics::abline (a = 0, b = 1, col = "blue")
        }
        else
        {
          dd = matrix (res [, , i], nrow = dim (res) [1])
          colnames (dd) = methodNames
          lim = c (min (dd), max (dd))
          graphics::pairs (dd, upper.panel = panel.compare, lower.panel = NULL, asp = 1, xlim = lim, ylim = lim)
        }
      }
    }
    else
      message ("Unknown evaluation")
  }

#' @keywords internal
# "Can. 1 (81.37 %)": the name of a canonical axis and the share of the trace it carries.
cda.axislabel <-
  function (x, axis)
  {
    eig = x$eig
    share = if (is.matrix (eig)) eig [axis, 2] else eig [2]
    name = colnames (x$proj) [axis]
    if (is.null (name))
      name = paste ("Can.", axis)
    return (paste0 (name, " (", round (share, 2), " %)"))
  }

#' Plot function for cda-class
#'
#' Plot the learning set (and test set) on the canonical axes obtained by Canonical Discriminant Analysis (function \code{CDA}).
#' @name plot.cda
#' @param x The classification model (object of class \code{cda-class}).
#' @param newdata The test set (\code{matrix} or \code{data.frame}).
#' @param axes The canonical axes to be printed (numeric \code{vector}). Ignored when there is
#' only one, which is what two classes give: the axis is then drawn against the observation
#' index, and the class means as horizontal lines.
#' @param legendpos Position of the legend, as in \code{\link{plotdata}}.
#' @param ... Other parameters, passed to the underlying plot.
#' @method plot cda
#' @export
#' @seealso \code{\link{CDA}}, \code{\link{predict.cda}}, \code{\link{cda-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = CDA (iris [, -5], iris [, 5])
#' plot (model)
plot.cda <-
  function (x, newdata = NULL, axes = 1:2, legendpos = "topleft", ...)
  {
    proj = as.matrix (x$proj)
    classes = factor (x$labels)
    levs = levels (classes)
    cols = 1 + seq_along (levs)
    flat = x$dim < 2
    test = NULL
    testclass = NULL
    if (!is.null (newdata))
    {
      test = as.matrix (cda.transform (x, newdata))
      testclass = factor (predict.cda (x, newdata), levels = levs)
    }
    if (flat)
    {
      # Two classes give a single canonical axis: there is nothing to plot it against but the
      # observation index, as plotdata() does for a dataset reduced to one variable. No aspect
      # ratio -- an index and a score share no unit.
      coord = cbind (Index = seq_len (nrow (proj)), proj [, 1])
      xlab = "Index"
      ylab = cda.axislabel (x, 1)
      asp = NA
      if (!is.null (test))
        test = cbind (Index = nrow (proj) + seq_len (nrow (test)), test [, 1])
    }
    else
    {
      axes = axes [1:2]
      coord = proj [, axes, drop = FALSE]
      xlab = cda.axislabel (x, axes [1])
      ylab = cda.axislabel (x, axes [2])
      asp = 1
      if (!is.null (test))
        test = test [, axes, drop = FALSE]
    }
    # One frame holding everything, so that the test set cannot fall outside it, then the
    # layers. The legend goes *inside* the plot, as in every other graphic of the package:
    # it used to sit in a band of its own, taking an eighth of the device whatever its size,
    # which truncated the entries on a wide screen.
    graphics::plot (rbind (coord, test), col = 0, xlab = xlab, ylab = ylab, asp = asp, ...)
    graphics::points (coord, col = cols [as.numeric (classes)], pch = 1)
    if (!is.null (test))
      graphics::points (test, col = cols [as.numeric (testclass)], pch = 3)
    if (flat)
      graphics::abline (h = tapply (proj [, 1], classes, mean), col = cols, lty = 2, lwd = 2)
    else
      graphics::points (t (sapply (levs, function (l)
                                   colMeans (coord [classes == l, , drop = FALSE]))),
                        col = cols, pch = 19, cex = 1.4)
    labels = c (levs, "class centre")
    lcol = c (cols, "grey40")
    lpch = c (rep (1, length (levs)), if (flat) NA else 19)
    # 0, not NA: legend() tests any (lty > 0), which NA turns into an error.
    llty = c (rep (0, length (levs)), if (flat) 2 else 0)
    if (!is.null (test))
    {
      labels = c (labels, "test set")
      lcol = c (lcol, "grey40")
      lpch = c (lpch, 3)
      llty = c (llty, 0)
    }
    graphics::legend (legendpos, legend = labels, col = lcol, pch = lpch, lty = llty,
                      bty = "n")
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
#' \dontrun{
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = BAGGING (d$train.x, d$train.y, NB)
#' predict (model, d$test.x)
#' model = ADABOOST (d$train.x, d$train.y, NB)
#' predict (model, d$test.x)
#' }
predict.boosting <- function (object, test, fuzzy = FALSE, ...)
{
  pred = lapply (object$models, function (model) predict (model, test, fuzzy, ...))
  weights = sapply (object$models, function (model) model$boostweight)
  res = NULL
  if (is.factor (object$y))
  {
    labels = levels (object$y)
    if (fuzzy)
    {
      # A *weighted* mean of the per-model probabilities, so that the weights AdaBoost spends
      # its run computing count here as they do in the hard vote below. BAGGING is unaffected,
      # all its weights being 1.
      total = sum (weights)
      res = Reduce ("+", Map (function (p, w) as.matrix (p) * w, pred, as.list (weights))) / total
      colnames (res) = labels
    }
    else
    {
      # The weighted vote, one pass per model over the whole test set rather than one
      # questionr::wtd.table() per observation -- which was 80% of the time spent predicting
      # with an ensemble (0.83 s for 100 models on 999 observations, against 0.02 s here).
      pred = sapply (pred, as.character)
      if (is.null (dim (pred)))
        pred = matrix (pred, nrow = 1)
      votes = matrix (0, nrow = nrow (pred), ncol = length (labels),
                      dimnames = list (NULL, labels))
      rows = seq_len (nrow (pred))
      for (j in seq_len (ncol (pred)))
      {
        cells = cbind (rows, match (pred [, j], labels))
        keep = !is.na (cells [, 2])
        votes [cells [keep, , drop = FALSE]] = votes [cells [keep, , drop = FALSE]] + weights [j]
      }
      res = factor (labels [max.col (votes, ties.method = "first")], levels = labels)
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

#' @keywords internal
predictmodel.cart <-
  function (object, test, fuzzy, ...)
  {
    if (is.vector (test))
    {
      test = matrix (test, ncol = 1)
      colnames (test) = "X"
    }
    if (fuzzy)
      res = stats::predict (object$model, as.data.frame (test))
    else
    {
      if (object$type == "class")
        res = stats::predict (object$model, as.data.frame (test), type = "class")
      else
        res = stats::predict (object$model, as.data.frame (test))
    }
    return (res)
  }

#' @keywords internal
predictmodel.lda <-
  function (object, test, fuzzy, ...)
  {
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    # drop = FALSE: see discriminant.variables(). With a single retained predictor,
    # test [, variables] would collapse to a bare vector, which MASS's predict methods reject
    # just as their fitting counterparts do.
    test = test [, object$variables, drop = FALSE]
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
    return (res)
  }

#' @keywords internal
predictmodel.lm <-
  function (object, test, fuzzy, ...)
  {
    if (is.vector (test))
      test = data.frame (X = test)
    return (stats::predict (object$model, test))
  }

#' @keywords internal
predictmodel.lr <-
  function (object, test, fuzzy, ...)
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
        # For a two-class problem, nnet::multinom returns the probability of the *second*
        # level only, as a plain vector: the first column of the result must therefore be
        # 1 - res, not res (same convention as predictmodel.xgb below).
        res = cbind (1 - res, res)
        colnames (res) = labels
      }
    }
    else
      res = stats::predict (object$model, test, type = "class")
    return (res)
  }

#' @keywords internal
predictmodel.mlp <-
  function (object, test, fuzzy, ...)
  {
    if (is.vector (test))
      test = data.frame (X = test)
    if (fuzzy)
    {
      res = stats::predict (object$model, test, type = "raw")
      labels = object$model$lev
      if (length (labels) == 2)
      {
        # Same as predictmodel.lr(): with two classes nnet::nnet has a single output unit,
        # holding the probability of the *second* level, so the first column is 1 - res.
        res = cbind (1 - res, res)
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
    return (res)
  }

#' @keywords internal
predictmodel.mlpreg <-
  function (object, test, fuzzy, ...)
  {
    rng = object$model$range [-1]
    rng [rng == 0] = 1 # Avoid division by zero for constant variables (NaN/Inf)
    d.norm = sweep (sweep (as.data.frame (test), 2, object$model$minimum [-1], FUN = "-"), 2, rng, FUN = "/")
    colnames (d.norm) = attr (attr (object$model$model$terms, "factor"), "dimnames") [[1]] [-1]
    return ((stats::predict (object$model$model, d.norm) * object$model$range [1]) + object$model$minimum [1])
  }

#' @keywords internal
predictmodel.mrv <-
  function (object, test, fuzzy, ...)
  {
    return (stats::predict (object$model, test, ncomp = object$model$ncomp))
  }

#' @keywords internal
predictmodel.nb <-
  function (object, test, fuzzy, ...)
  {
    type = "class"
    if (fuzzy)
      type = "raw"
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    return (stats::predict (object$model, test, type = type))
  }

#' @keywords internal
predictmodel.regularization <-
  function (object, test, fuzzy, ...)
  {
    return (stats::predict (object$model, as.matrix (test)))
  }

#' @keywords internal
predictmodel.rf <-
  function (object, test, fuzzy, ...)
  {
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    # randomForest() needs type = "prob" to return class-membership scores; without it the
    # fuzzy = TRUE request was silently ignored and hard labels came back, which left
    # performance (type = "roc") with no score to draw from.
    if (fuzzy)
      return (stats::predict (object$model, test, type = "prob"))
    return (stats::predict (object$model, test))
  }

#' @keywords internal
predictmodel.svm <-
  function (object, test, fuzzy, ...)
  {
    if (fuzzy)
    {
      res = attr (stats::predict (object$model, test, probability = TRUE), "probabilities")
      # e1071 orders those columns by its own internal class numbering, not by the levels of
      # the target. Every other method of the package returns them in level order, and
      # predict (model, x, fuzzy = TRUE) [, 2] must mean the same class whichever produced it.
      lev = object$model$levels
      if ((!is.null (lev)) && all (lev %in% colnames (res)))
        res = res [, lev, drop = FALSE]
    }
    else
      res = stats::predict (object$model, test, probability = FALSE)
    return (res)
  }

#' @keywords internal
predictmodel.xgb <-
  function (object, test, fuzzy, ...)
  {
    if (fuzzy)
    {
      # type = "response" returns a matrix of per-class probabilities for multi-class
      # objectives, but a plain vector (probability of the *last* factor level) for binary
      # ones; the latter is reshaped into the 2-column matrix used everywhere else in the
      # package for fuzzy/probabilistic predictions.
      res = stats::predict (object$model, as.matrix (test), type = "response")
      if (!is.matrix (res))
        res = cbind (1 - res, res)
      colnames (res) = object$lev
    }
    else
      # type = "class" already returns a factor using the levels seen during training.
      res = stats::predict (object$model, as.matrix (test), type = "class")
    return (res)
  }


#' Regression using Gradient Boosting
#'
#' This function builds a regression model using Gradient Boosting. It is the regression
#' counterpart of \code{\link{GRADIENTBOOSTING}}, which classifies.
#' @name GBREG
#' @param x Predictor values of the training set, as a \code{matrix} or \code{data.frame}.
#' @param y Target values of the training set (a numeric \code{vector}).
#' @param ntree The number of trees in the ensemble.
#' @param learningrate The learning rate (between 0 and 1).
#' @param seed A specified seed for random number generation (row/column subsampling, if used
#' via \code{...}).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{GBREG} does not yet
#' implement hyperparameter tuning.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{GBREG} does not produce a plot.
#' @param ... Other parameters, passed to \code{\link[xgboost]{xgboost}}.
#' @return The regression model.
#' @export
#' @seealso \code{\link{GRADIENTBOOSTING}}, \code{\link{LINREG}}, \code{\link{SVR}},
#' \code{\link[xgboost]{xgboost}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (trees)
#' d = splitdata (trees, 3)
#' model = GBREG (d$train.x, d$train.y)
#' evaluation (predict (model, d$test.x), d$test.y)
#' }
GBREG <-
  function (x, y, ntree = 500, learningrate = 0.3,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.factor (y))
        stop ("GBREG: the target must be numeric -- this is the regression counterpart of ",
              "GRADIENTBOOSTING(), which is the one to use on class labels.")
      setseed (seed)
      # As in GRADIENTBOOSTING(): since xgboost >= 2.1 the objective is inferred from 'y', a
      # numeric one selecting the squared-error regression objective.
      model = xgboost::xgboost (x = as.matrix (x), y = as.numeric (y), nrounds = ntree,
                                learning_rate = learningrate, verbosity = 0, ...)
      res = list (model = model, method = "XGBREG")
      class (res) = "model"
    }
    return (res)
  }

#' @keywords internal
predictmodel.xgbreg <-
  function (object, test, fuzzy, ...)
    as.vector (stats::predict (object$model, as.matrix (test)))

#' @keywords internal
predictmodel.penalizedlr <-
  function (object, test, fuzzy, ...)
  {
    if (fuzzy)
    {
      p = stats::predict (object$model, as.matrix (test), type = "response")
      # glmnet returns an n x K x 1 array for a multinomial family, and an n x 1 matrix (the
      # probability of the *second* level) for a binomial one.
      if (length (dim (p)) == 3)
        p = p [, , 1]
      else
        p = cbind (1 - as.vector (p), as.vector (p))
      colnames (p) = object$lev
      return (p)
    }
    return (factor (as.vector (stats::predict (object$model, as.matrix (test), type = "class")),
                    levels = object$lev))
  }

#' @keywords internal
# Lookup table for the predictmodel.* dispatch used by predict.model(). Built lazily (as a
# function rather than a static list) so it does not depend on file collation order. LDA and
# QDA share the same implementation, as in the original if/else chain.
predictmodel.functions <-
  function ()
  {
    list (CART = predictmodel.cart,
          LDA = predictmodel.lda,
          QDA = predictmodel.lda,
          lm = predictmodel.lm,
          LR = predictmodel.lr,
          MLP = predictmodel.mlp,
          MLPREG = predictmodel.mlpreg,
          MRV = predictmodel.mrv,
          NB = predictmodel.nb,
          regularization = predictmodel.regularization,
          RF = predictmodel.rf,
          penalizedlr = predictmodel.penalizedlr,
          SVM = predictmodel.svm,
          XGB = predictmodel.xgb,
          XGBREG = predictmodel.xgbreg)
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
    funs = predictmodel.functions ()
    fun = funs [[object$method]]
    if (is.null (fun))
      res = stats::predict (object$model, test, ...)
    else
      res = fun (object, test, fuzzy, ...)
    return (res)
  }

#' @keywords internal
# Shared fit + predict logic for a single train/test split, used by every protocol.*
# function below (bootstrap, crossvalidation, holdout, loocv, train). Handles both the
# single-method and multi-method (list of methods) cases uniformly, so the split-generation
# logic specific to each protocol is the only thing left duplicated.
protocol.fitpredict <-
  function (methods, indices, learn, learny, test, methodparameters, fuzzy = FALSE, ...)
  {
    # Fitting one model per split is no place for a plot, so 'graph = FALSE' is imposed --
    # unless the caller asked for something else. Passing it unconditionally meant that
    # performance (LINREG, ..., graph = FALSE), i.e. the obvious way to silence a method that
    # draws by default, stopped on 'formal argument "graph" matched by multiple actual
    # arguments'.
    extra = list (...)
    if (is.null (extra$graph))
      extra$graph = FALSE
    fit = function (method, params)
      do.call (method, c (list (learn, learny), list (methodparameters = params), extra))
    if (length (methods) == 1)
      models = list (fit (methods, methodparameters))
    else
      models = lapply (indices, function (i) fit (methods [[i]], methodparameters [[i]]))
    if (fuzzy)
      # One probability matrix per method, always as a list -- performance() stacks the
      # matrices of the successive splits before extracting the positive class column.
      return (lapply (models, function (model) as.matrix (stats::predict (model, test, fuzzy = TRUE, ...))))
    if (length (models) == 1)
      return (stats::predict (models [[1]], test, ...))
    return (sapply (models, function (model) as.numeric (stats::predict (model, test, ...))))
  }

#' @keywords internal
# Every protocol.* function below takes the *same* argument list, whether it uses all of it or
# not, because performance() passes all of it by name. An argument left to fall into '...'
# would travel on to the learning method, and into randomForest(), nnet(), svm() or xgboost(),
# which swallow unknown arguments silently.
protocol.arguments <-
  c ("methods", "train.x", "train.y", "test.x", "test.y", "train.size",
     "methodparameters", "nruns", "nfolds", "seed", "fuzzy")

#' @keywords internal
protocol.bootstrap <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL, train.size = NULL,
            methodparameters = NULL, nruns = 10, nfolds = 10, seed = NULL, fuzzy = FALSE,
            stratify = TRUE, ...)
  {
    setseed (seed)
    predictions = NULL
    targets = NULL
    n = length (train.y)
    indices = 1:length (methods)
    samples = matrix (sample (n, n * nruns, replace = TRUE), ncol = nruns)
    for (i in 1:nruns)
    {
      s = samples [, i]
      targets = c (targets, list (train.y [-s]))
      learn = train.x [s, ]
      if (is.vector (learn))
        learn = data.frame (X = learn)
      test = train.x [-s, ]
      if (is.vector (test))
        test = data.frame (X = test)
      rownames (learn) = 1:nrow (learn)
      pred = protocol.fitpredict (methods, indices, learn, train.y [s], test, methodparameters, fuzzy = fuzzy, ...)
      predictions = c (predictions, list (pred))
    }
    return (list (predictions = predictions, targets = targets))
  }

#' @keywords internal
protocol.crossvalidation <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL, train.size = NULL,
            methodparameters = NULL, nruns = 10, nfolds = 10, seed = NULL, fuzzy = FALSE,
            stratify = TRUE, ...)
  {
    setseed (seed)
    predictions = NULL
    targets = NULL
    n = length (train.y)
    indices = 1:length (methods)
    for (i in 1:nruns)
    {
      # 'seed + i' is numeric (0) when seed is NULL; setseed() treats that as "no seed given".
      setseed (seed + i)
      # Stratified folds keep the class proportions of the whole sample in every fold. Without
      # them a rare class can be missing from a fold entirely -- so the model is asked to
      # predict a class it never saw, or scored on a fold that contains none of it.
      if (stratify && stratifiable (train.y))
      {
        s = seq_len (n)
        folds = stratified.folds (train.y, nfolds)
      }
      else
      {
        s = sample (n, n, replace = FALSE)
        folds = rep (1:nfolds, diff (round ((n / nfolds) * 0:(nfolds))))
      }
      for (j in 1:nfolds)
      {
        slearn = s [folds != j]
        stest = s [folds == j]
        targets = c (targets, list (train.y [stest]))
        learn = train.x [slearn, ]
        if (is.vector (learn))
          learn = data.frame (X = learn)
        test = train.x [stest, ]
        if (is.vector (test))
          test = data.frame (X = test)
        rownames (learn) = 1:nrow (learn)
        pred = protocol.fitpredict (methods, indices, learn, train.y [slearn], test, methodparameters, fuzzy = fuzzy, ...)
        predictions = c (predictions, list (pred))
      }
    }
    return (list (predictions = predictions, targets = targets))
  }

#' @keywords internal
protocol.holdout <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL,
            train.size = round (0.7 * length (train.y)),
            methodparameters = NULL, nruns = 10, nfolds = 10, seed = NULL, fuzzy = FALSE,
            stratify = TRUE, ...)
  {
    setseed (seed)
    if (is.null (test.x) | is.null (test.y))
    {
      if (train.size < 1)
        train.size = round (train.size * length (train.y))
      if (stratify && stratifiable (train.y))
        s = stratified.sample (train.y, train.size)
      else
        s = sample (nrow (train.x), train.size)
      test.x = train.x [-s, ]
      train.x = train.x [s, ]
      if (is.vector (train.x))
        train.x = data.frame (X = train.x)
      if (is.vector (test.x))
        test.x = data.frame (X = test.x)
      test.y = train.y [-s]
      train.y = train.y [s]
    }
    setseed (seed)
    indices = 1:length (methods)
    predictions = protocol.fitpredict (methods, indices, train.x, train.y, test.x, methodparameters, fuzzy = fuzzy, ...)
    targets = test.y
    return (list (predictions = list (predictions), targets = list (targets)))
  }

#' @keywords internal
protocol.loocv <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL, train.size = NULL,
            methodparameters = NULL, nruns = 10, nfolds = 10, seed = NULL, fuzzy = FALSE,
            stratify = TRUE, ...)
  {
    setseed (seed)
    targets = train.y
    indices = 1:length (methods)
    fit = function (j)
      protocol.fitpredict (methods, indices, train.x [-j, , drop = FALSE], train.y [-j],
                           train.x [j, , drop = FALSE], methodparameters, fuzzy = fuzzy, ...)
    if (fuzzy)
    {
      # fit() returns one probability matrix (of a single row here) per method; stack them
      # method by method so the result has the same shape as the other protocols'.
      byobs = lapply (1:length (train.y), fit)
      predictions = lapply (seq_len (length (methods)),
                            function (m) do.call (rbind, lapply (byobs, function (obs) obs [[m]])))
    }
    else if (length (methods) == 1)
      predictions = sapply (1:length (train.y), fit)
    else
      predictions = t (sapply (1:length (train.y), fit))
    return (list (predictions = list (predictions), targets = list (targets)))
  }

#' @keywords internal
protocol.train <-
  function (methods, train.x, train.y, test.x = NULL, test.y = NULL, train.size = NULL,
            methodparameters = NULL, nruns = 10, nfolds = 10, seed = NULL, fuzzy = FALSE,
            stratify = TRUE, ...)
  {
    setseed (seed)
    indices = 1:length (methods)
    predictions = protocol.fitpredict (methods, indices, train.x, train.y, train.x, methodparameters, fuzzy = fuzzy, ...)
    targets = train.y
    return (list (predictions = list (predictions), targets = list (targets)))
  }

#' @keywords internal
# Lookup table for the protocol.* dispatch used by performance(). Built lazily (as a function
# rather than a static list) so it does not depend on file collation order. This was the last
# get (paste (...)) based dispatch left in the package.
protocol.functions <-
  function ()
  {
    list (bootstrap = protocol.bootstrap,
          crossvalidation = protocol.crossvalidation,
          holdout = protocol.holdout,
          loocv = protocol.loocv,
          train = protocol.train)
  }

#' Classification using Quadratic Discriminant Analysis
#'
#' This function builds a classification model using Quadratic Discriminant Analysis.
#' @name QDA
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{QDA} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{QDA} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link[MASS]{qda}}
#' @examples
#' require (datasets)
#' data (iris)
#' QDA (iris [, -5], iris [, 5])
QDA <-
  function (train, labels, tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      variables = discriminant.variables (train, "QDA")
      model = MASS::qda (x = train [, variables, drop = FALSE], grouping = labels)
      res = list (model = model, method = "QDA", variables = variables)
      class (res) = "model"
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
#' @inheritParams tune.doc
#' @param seed A specified seed for random number generation (bootstrap sampling of the trees
#' and, when \code{nvar} is smaller than the total number of variables, the candidate variables
#' drawn at each split).
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{RANDOMFOREST} does
#' not yet implement hyperparameter tuning, so \code{tune = TRUE} returns an empty
#' \code{params} object and there is nothing for \code{methodparameters} to override.
#' @param graph Present for interface consistency with \code{\link{performance}} (which
#' always passes it when fitting a model). Currently unused: \code{RANDOMFOREST} does not
#' produce a plot.
#' @param ... Other parameters, forwarded to \code{\link[randomForest]{randomForest}}.
#' @return The classification model.
#' @export
#' @seealso \code{\link[randomForest]{randomForest}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' RANDOMFOREST (iris [, -5], iris [, 5])
#' }
RANDOMFOREST <-
  function (train, labels,
            ntree = 500,
            nvar = if (!is.null (labels) && !is.factor (labels)) max (floor (ncol (train)/3), 1) else floor (sqrt (ncol (train))),
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      setseed (seed)
      res = randomForest::randomForest(x = train, y = labels, ntree = ntree, mtry = nvar, ...)
      res = list (model = res, method = "RF")
      class (res) = "model"
    }
    return (res)
  }

#' Plot ROC Curves
#'
#' This function plots ROC Curves of one or several classification predictions.
#'
#' A ROC curve needs a \emph{score}: the higher it is, the more likely the observation is to
#' belong to the positive class. The natural one is the estimated probability returned by
#' \code{predict (model, x, fuzzy = TRUE)}. Hard class labels give only two distinct values,
#' so the "curve" reduces to three points and the area under it says very little. That coarse
#' version is available on purpose -- it makes a useful comparison -- but it has to be asked
#' for: pass hard labels, or \code{type = "hard"}.
#' @name roc.curves
#' @param predictions The predictions of one or several classification models. Four shapes are
#' accepted: a \code{factor} of hard labels (one model); a numeric \code{vector} of scores for
#' the positive class (one model); a \code{matrix} of class probabilities, i.e. one column per
#' class named after it, as returned by \code{predict (model, x, fuzzy = TRUE)} (one model);
#' or any other \code{matrix}/\code{data.frame}, read as one column per model.
#' @param gt Actual labels of the dataset (\code{factor} or \code{vector}), two classes only.
#' @param methods.names The name of the compared methods (\code{vector}).
#' @param positive The label of the positive class. Defaults to the first level of \code{gt},
#' as everywhere else in the package -- note that for the usual alphabetical ordering this is
#' often the \emph{negative} class, so it is worth setting explicitly.
#' @param type \code{"auto"} (default) reads \code{predictions} according to its shape,
#' \code{"fuzzy"} requires scores and refuses hard labels, \code{"hard"} reduces everything to
#' the predicted class first (this is the coarse three-point curve discussed above).
#' @param ... Other parameters, passed to the underlying plot.
#' @return Nothing; the curves are drawn on the current graphics device.
#' @export
#' @seealso \code{\link{cost.curves}}, \code{\link{performance}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = iris
#' levels (d [, 5]) = c ("+", "+", "-") # Building a two classes dataset
#' model.nb = NB (d [, -5], d [, 5])
#' model.lda = LDA (d [, -5], d [, 5])
#' # From the estimated probabilities: the meaningful curve
#' roc.curves (predict (model.nb, d [, -5], fuzzy = TRUE), d [, 5], positive = "+")
#' # Two models compared, one score column each
#' roc.curves (cbind (NB = predict (model.nb, d [, -5], fuzzy = TRUE) [, "+"],
#'                    LDA = predict (model.lda, d [, -5], fuzzy = TRUE) [, "+"]),
#'             d [, 5], c ("NB", "LDA"), positive = "+")
#' # The same predictions reduced to hard labels: three points, and little to read
#' roc.curves (predict (model.nb, d [, -5]), d [, 5], positive = "+", type = "hard")
roc.curves <-
  function (predictions, gt, methods.names = NULL, positive = levels (factor (gt)) [1],
            type = c ("auto", "fuzzy", "hard"), ...)
  {
    curve.plot (predictions, gt, methods.names, positive, type,
                measure = c ("tpr", "fpr"), asp = 1, ...)
  }

#' Classification using one-level decision tree
#'
#' This function builds a classification model using CART with maxdepth = 1.
#' @name STUMP
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param randomvar If \code{TRUE}, the stump is built on a single variable drawn at random
#' instead of the best one (useful to build weak learners for an ensemble method). Note that
#' the model then differs from one call to the next unless \code{seed} is set. Defaults to
#' \code{FALSE}, i.e. the usual decision stump, split on the variable selected by CART.
#' @inheritParams tune.doc
#' @param seed A specified seed for random number generation (used only if \code{randomvar} is \code{TRUE}).
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{STUMP} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{STUMP} does not produce a plot.
#' @param ... Other parameters.
#' @return The classification model.
#' @export
#' @seealso \code{\link{CART}}
#' @examples
#' require (datasets)
#' data (iris)
#' STUMP (iris [, -5], iris [, 5])
#' STUMP (iris [, -5], iris [, 5], randomvar = TRUE, seed = 0)
STUMP <-
  function (train, labels, randomvar = FALSE,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    new = train
    if (randomvar && (!is.vector (train)))
    {
      setseed (seed)
      var = sample (ncol (train), 1)
      new = matrix (train [, var], ncol = 1)
      colnames (new) = colnames (train) [var]
    }
    return (CART (new, labels, minsplit = 1, maxdepth = 1, cp = 0, tune = tune, ...))
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
#' @inheritParams tune.doc
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{SVM} does not produce a plot.
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVMl}}, \code{\link{SVMr}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' SVM (iris [, -5], iris [, 5], kernel = "linear", cost = 1)
#' SVM (iris [, -5], iris [, 5], kernel = "radial", gamma = 1, cost = 1)
#' }
SVM <-
  function (train,
            labels,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            kernel = c ("radial", "linear"),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    check.numeric.predictors (train, "SVM")
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
    }
    if (kernel [1] == "linear")
      gamma = 0
    if (length (gamma) > 1 | length (cost) > 1)
    {
      # No 'probability = TRUE' while searching the grid: e1071 fits the probability model
      # (Platt scaling) by an *internal* five-fold cross-validation, which turns every fit of
      # the grid into six, all discarded with the grid point they belong to. Only the retained
      # model is ever asked for probabilities.
      tuned = e1071::tune.svm (train, labels,
                               gamma = gamma, cost = cost, kernel = kernel [1],
                               tunecontrol = tune.scheme (nfolds), ...)
      # Always refitted, rather than reading tune.svm()'s own 'best.model': that one carries no
      # probability model, and is empty altogether whenever tune.svm() could not refit it.
      model = e1071::svm (train, labels, gamma = tuned$best.parameters$gamma,
                          cost = tuned$best.parameters$cost, kernel = kernel [1],
                          probability = TRUE, ...)
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
#' @inheritParams tune.doc
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVM}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' SVMl (iris [, -5], iris [, 5], cost = 1)
#' }
SVMl <-
  function (train,
            labels,
            cost = 2^(-3:3),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    return (SVM (
      train = train,
      labels = labels,
      gamma = NULL,
      cost = cost,
      kernel = "linear",
      nfolds = nfolds,
      tune = tune,
      methodparameters = methodparameters,
      graph = graph,
      seed = seed,
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
#' @inheritParams tune.doc
#' @param ... Other arguments.
#' @return The classification model.
#' @export
#' @seealso \code{\link[e1071]{svm}}, \code{\link{SVM}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' SVMr (iris [, -5], iris [, 5], gamma = 1, cost = 1)
#' }
SVMr <-
  function (train,
            labels,
            gamma = 2^(-3:3),
            cost = 2^(-3:3),
            nfolds = 10,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    return (SVM (
      train = train,
      labels = labels,
      gamma = gamma,
      cost = cost,
      kernel = "radial",
      nfolds = nfolds,
      tune = tune,
      methodparameters = methodparameters,
      graph = graph,
      seed = seed,
      ...
    ))
  }

#' @keywords internal
# Shared layout of the print methods below: a title line, then aligned "field: value" rows.
# Keeps every method to a description of what it holds rather than to formatting code.
shortprint <-
  function (title, fields)
  {
    cat (title, "\n", sep = "")
    fields = fields [!sapply (fields, is.null)]
    if (length (fields) == 0)
      return (invisible (NULL))
    width = max (nchar (names (fields)))
    for (n in names (fields))
      cat ("  ", formatC (n, width = width, flag = "-"), " : ",
           paste (fields [[n]], collapse = ", "), "\n", sep = "")
    invisible (NULL)
  }

#' @keywords internal
# "150 observations, 4 variables", or "150 observations" for a bare vector.
sizelabel <-
  function (x)
  {
    if (is.null (x))
      return (NULL)
    if (is.null (dim (x)))
      return (paste (length (x), "observations"))
    return (paste (nrow (x), "observations,", ncol (x), "variables"))
  }

#' @keywords internal
# "3 classes (setosa, versicolor, virginica)" or, for a numeric target, its range.
targetlabel <-
  function (y)
  {
    if (is.null (y))
      return (NULL)
    if (is.factor (y) || is.character (y))
    {
      y = factor (y)
      return (paste0 (nlevels (y), " classes (", paste (levels (y), collapse = ", "), ")"))
    }
    return (paste0 ("numeric, from ", signif (min (y), 4), " to ", signif (max (y), 4)))
  }

#' Print a classification or regression model
#'
#' Prints a short description of the model -- the method it was obtained with, and the size of
#' the data it was fitted on -- instead of dumping the underlying list. Use
#' \code{summary (model)} for the full detail of the wrapped model.
#' @name print.model
#' @param x The model to be printed (object of class \code{\link{model-class}}).
#' @param ... Other parameters.
#' @export
#' @method print model
#' @seealso \code{\link{model-class}}, \code{\link{summary.model}}, \code{\link{predict.model}}
#' @examples
#' require (datasets)
#' data (iris)
#' NB (iris [, -5], iris [, 5])
print.model <-
  function (x, ...)
  {
    shortprint (paste0 ("Classification/regression model (", x$method, ")"),
                  list ("method" = x$method,
                        "wrapped model" = paste (class (x$model), collapse = ", ")))
  }

#' Summary of a classification or regression model
#'
#' Prints the short description of \code{\link{print.model}}, followed by the summary of the
#' wrapped model itself.
#' @name summary.model
#' @param object The model (object of class \code{\link{model-class}}).
#' @param ... Other parameters, passed to the summary of the wrapped model.
#' @export
#' @method summary model
#' @seealso \code{\link{model-class}}, \code{\link{print.model}}
#' @examples
#' require (datasets)
#' data (iris)
#' summary (NB (iris [, -5], iris [, 5]))
summary.model <-
  function (object, ...)
  {
    print (object)
    cat ("\n")
    return (summary (object$model, ...))
  }

#' Print a training/test split
#'
#' Prints the sizes of the two parts of a split and the target they share.
#' @name print.dataset
#' @param x The split (object of class \code{\link{dataset-class}}, created by
#' \code{\link{splitdata}}).
#' @param ... Other parameters.
#' @export
#' @method print dataset
#' @seealso \code{\link{dataset-class}}, \code{\link{splitdata}}
#' @examples
#' require (datasets)
#' data (iris)
#' splitdata (iris, 5)
print.dataset <-
  function (x, ...)
  {
    n = length (x$train.y) + length (x$test.y)
    shortprint ("Training/test split",
                  list ("training set" = paste0 (sizelabel (x$train.x), " (",
                                                 round (100 * length (x$train.y) / n), "%)"),
                        "test set" = paste0 (sizelabel (x$test.x), " (",
                                             round (100 * length (x$test.y) / n), "%)"),
                        "target" = targetlabel (x$train.y),
                        "class balance" = if (is.factor (x$train.y))
                          paste0 ("train ", paste (as.vector (table (x$train.y)), collapse = "/"),
                                  ", test ", paste (as.vector (table (x$test.y)), collapse = "/"))
                        else NULL))
  }

#' Print a feature selection result
#'
#' Prints which features were selected, by which algorithm and criteria, instead of dumping the
#' underlying list.
#' @name print.selection
#' @param x The result (object of class \code{\link{selection-class}}, created by
#' \code{\link{selectfeatures}}).
#' @param ... Other parameters.
#' @export
#' @method print selection
#' @seealso \code{\link{selection-class}}, \code{\link{selectfeatures}}
#' @examples
#' require (datasets)
#' data (iris)
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "cfs")
print.selection <-
  function (x, ...)
  {
    selected = if (!is.null (x$features)) x$features else x$selection
    shortprint ("Feature selection",
                  list ("algorithm" = x$algorithm,
                        "univariate criterion" = x$univariate,
                        "multivariate criterion" = x$multivariate,
                        "features kept" = paste0 (length (x$selection), ": ",
                                                  paste (selected, collapse = ", ")),
                        "score" = if (!is.null (x$multieval)) signif (x$multieval, 4) else NULL))
  }

#' Print tuned method parameters
#'
#' Prints the hyperparameters a method retained, or says that it has none.
#' @name print.params
#' @param x The parameters (object of class \code{\link{params-class}}, obtained by calling a
#' classification method with \code{tune = TRUE}).
#' @param ... Other parameters.
#' @export
#' @method print params
#' @seealso \code{\link{params-class}}, \code{\link{performance}}
#' @examples
#' require (datasets)
#' data (iris)
#' # A small grid, so that the example stays fast; the defaults search a much larger one.
#' SVM (iris [, -5], iris [, 5], gamma = 2^(-2:0), cost = 2^(0:2), tune = TRUE)
#' # A method with nothing to tune says so.
#' NB (iris [, -5], iris [, 5], tune = TRUE)
print.params <-
  function (x, ...)
  {
    if (length (x) == 0)
      return (shortprint ("Method parameters", list ("none" = "this method has none to tune")))
    values = lapply (unclass (x), function (v) paste (signifnum (v), collapse = ", "))
    shortprint ("Method parameters", values)
  }

#' @keywords internal
signifnum <-
  function (v) if (is.numeric (v)) signif (v, 4) else v

#' Print a canonical discriminant analysis
#'
#' Prints the size of the analysis and the variance carried by its axes.
#' @name print.cda
#' @param x The model (object of class \code{\link{cda-class}}, created by \code{\link{CDA}}).
#' @param ... Other parameters.
#' @export
#' @method print cda
#' @seealso \code{\link{cda-class}}, \code{\link{CDA}}, \code{\link{plot.cda}}
#' @examples
#' require (datasets)
#' data (iris)
#' CDA (iris [, -5], iris [, 5])
print.cda <-
  function (x, ...)
  {
    eig = x$eig
    variance = if (is.matrix (eig)) eig [, 2] else eig [2]
    power = if (is.matrix (eig)) eig [, 4] else eig [4]
    shortprint ("Canonical discriminant analysis",
                  list ("training set" = sizelabel (x$train),
                        "target" = targetlabel (x$labels),
                        "canonical axes" = x$dim,
                        "variance (%)" = paste (round (variance, 2), collapse = ", "),
                        "discriminant power (%)" = paste (round (power, 2), collapse = ", ")))
  }

#' Print an ensemble model
#'
#' Prints how many models the ensemble holds and what it was fitted on.
#' @name print.boosting
#' @param x The model (object of class \code{\link{boosting-class}}, created by
#' \code{\link{ADABOOST}} or \code{\link{BAGGING}}).
#' @param ... Other parameters.
#' @export
#' @method print boosting
#' @seealso \code{\link{boosting-class}}, \code{\link{ADABOOST}}, \code{\link{BAGGING}}
#' @examples
#' require (datasets)
#' data (iris)
#' BAGGING (iris [, -5], iris [, 5], LDA, nsamples = 5)
print.boosting <-
  function (x, ...)
  {
    kept = length (x$models)
    shortprint ("Ensemble model",
                  list ("models" = if ((!is.null (x$nsamples)) && (kept < x$nsamples))
                                     paste0 (kept, " (of the ", x$nsamples, " asked for: ",
                                             "boosting stopped early)")
                                   else kept,
                        "training set" = sizelabel (x$x),
                        "target" = targetlabel (x$y)))
  }

#' Print a \emph{K}-nearest-neighbours model
#'
#' Prints the size of the training set the model memorised, and the number of neighbours used.
#' @name print.knn
#' @param x The model (object of class \code{\link{knn-class}}, created by \code{\link{KNN}}).
#' @param ... Other parameters.
#' @export
#' @method print knn
#' @seealso \code{\link{knn-class}}, \code{\link{KNN}}
#' @examples
#' require (datasets)
#' data (iris)
#' KNN (iris [, -5], iris [, 5])
print.knn <-
  function (x, ...)
  {
    shortprint ("K-nearest-neighbours model",
                  list ("neighbours (k)" = x$k,
                        "training set" = sizelabel (x$train),
                        "target" = targetlabel (x$labels)))
  }
