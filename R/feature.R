#' Feature selection
#'
#' This class contains the result of feature selection algorithms.
#' @name selection-class
#' @slot selection A vector of integers indicating the selected features.
#' @slot unieval The evaluation of the features (univariate).
#' @slot multieval The evaluation of the selected features (multivariate).
#' @slot algorithm The algorithm used to select features.
#' @slot univariate The evaluation criterion (univariate).
#' @slot nbfeatures The number of features to be kept.
#' @slot threshold The threshold to decide whether a feature is kept or not..
#' @slot multivariate The evaluation criterion (multivariate).
#' @slot dataset The dataset described by the selected features only.
#' @slot model The classification model.
#' @exportClass selection
#' @seealso \code{\link{FEATURESELECTION}}, \code{\link{predict.selection}}, \code{\link{selectfeatures}}
setClass ("selection",
          representation (selection = "vector",
                          unieval = "vector",
                          multieval = "numeric",
                          algorithm = "character",
                          univariate = "character",
                          nbfeatures = "numeric",
                          threshold = "numeric",
                          multivariate = "character",
                          dataset = "matrix",
                          model = "ANY"))

#' Classification with Feature selection
#'
#' Apply a classification method after a subset of features has been selected.
#' @name FEATURESELECTION
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param algorithm The feature selection algorithm.
#' @param unieval The (univariate) evaluation criterion. \code{uninb}, \code{unithreshold} or \code{multieval} must be specified.
#' @param uninb The number of selected feature (univariate evaluation).
#' @param unithreshold The threshold for selecting feature (univariate evaluation).
#' @param multieval The (multivariate) evaluation criterion.
#' @param wrapmethod The classification method used for the wrapper evaluation.
#' @param mainmethod The final method used for data classification. If a wrapper evaluation is used, the same classification method should be used.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{selectfeatures}}, \code{\link{predict.selection}}, \code{\link{selection-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' FEATURESELECTION (iris [, -5], iris [, 5], uninb = 2, mainmethod = LDA)
FEATURESELECTION <-
  function (train,
            labels,
            algorithm = c ("ranking", "forward", "backward", "exhaustive"),
            unieval = if (algorithm [1] == "ranking") c ("fisher", "fstat", "relief", "inertiaratio") else NULL,
            uninb = NULL,
            unithreshold = NULL,
            multieval = if (algorithm [1] == "ranking") NULL else c ("cfs", "fstat", "inertiaratio", "wrapper"),
            wrapmethod = NULL,
            mainmethod = wrapmethod,
            ...)
  {
    selection = selectfeatures (train, labels, algorithm, unieval, uninb, unithreshold, multieval, wrapmethod, keep = TRUE, ...)
    selection$model = mainmethod (selection$dataset, labels, ...)
    return (selection)
  }

#' @keywords internal
fs.backward <-
  function (train, labels, multieval, ...)
  {
    select = rep (TRUE, ncol (train))
    besteval = multieval (train, labels, type = "multivariate", ...)
    bestselect = select
    while (sum (select) > 1)
    {
      removenext = which (select)
      eval = sapply (removenext, function (index)
      {
        nextselect = select
        nextselect [index] = FALSE
        subset = train [, nextselect]
        return (multieval (subset, labels, type = "multivariate", ...))
      })
      localbesteval = max (eval)
      select [removenext [which.max (eval)]] = FALSE
      if (localbesteval > besteval)
      {
        besteval = localbesteval
        bestselect = select
      }
    }
    res = list (selection = which (bestselect), multieval = max (besteval), algorithm = "backward")
    return (res)
  }

#' @keywords internal
fs.exhaustive <-
  function (train, labels, multieval, ...)
  {
    indices = 1:(2^ncol (train) - 1)
    eval = sapply (indices, function (index)
    {
      select = as.logical (as.integer (intToBits (index) [1:ncol (train)]))
      subset = train [, select]
      return (multieval (subset, labels, type = "multivariate", ...))
    })
    res = as.logical (as.integer (intToBits (which.max (eval)) [1:ncol (train)]))
    res = list (selection = which (res), multieval = max (eval), algorithm = "exhaustive")
    return (res)
  }

#' @keywords internal
fs.forward <-
  function (train, labels, multieval, ...)
  {
    select = rep (FALSE, ncol (train))
    besteval = -Inf
    bestselect = select
    while (any (!select))
    {
      addnext = which (!select)
      eval = sapply (addnext, function (index)
      {
        nextselect = select
        nextselect [index] = TRUE
        subset = train [, nextselect]
        return (multieval (subset, labels, type = "multivariate", ...))
      })
      localbesteval = max (eval)
      select [addnext [which.max (eval)]] = TRUE
      if (localbesteval > besteval)
      {
        besteval = localbesteval
        bestselect = select
      }
    }
    res = list (selection = which (bestselect), multieval = max (besteval), algorithm = "forward")
    return (res)
  }

#' @keywords internal
fs.ranking <-
  function (train, labels, unieval, uninb, unithreshold, multieval, wrapmethod, ...)
  {
    eval = unieval (train, labels, type = "univariate")
    res = NULL
    if (!is.null (uninb) && (uninb > 0))
    {
      selection = order (eval, decreasing = TRUE) [1:uninb]
      res = list (selection = selection, unieval = eval, algorithm = "ranking", nbfeatures = uninb)
    }
    else if (!is.null (unithreshold) && (unithreshold < max (eval)))
    {
      selection = which (eval > unithreshold)
      names (selection) = NULL
      res = list (selection = selection, unieval = eval, algorithm = "ranking", threshold = unithreshold)
    }
    else if ((!is.null (multieval)) && (!is.null (wrapmethod)))
    {
      size = 1:ncol (train)
      features = order (eval, decreasing = TRUE)
      meval = sapply (size, function (index)
      {
        subset = train [, features [1:index]]
        return (multieval (subset, labels, type = "multivariate", wrapmethod = wrapmethod, ...))
      })
      selection = sort (features [1:which.max (meval)])
      res = list (selection = selection, unieval = eval, multieval = max (meval), algorithm = "ranking")
    }
    else
      message ("Cannot select features")
    return (res)
  }

#' @keywords internal
fseval.cfs <-
  function (train, labels, type = c ("multivariate", "univariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (type [1] == "univariate")
    {
      message ("CFS is an multivariate measure")
      return (NULL)
    }
    else
    {
      k = ncol (train)
      r = 1
      if (k > 1)
      {
        r = stats::cor (train)
        r = mean (abs (r [lower.tri(r)]))
      }
      return (k * mean (fseval.relief (train, labels)) / sqrt (k + (k * (k - 1) * r)))
    }
  }

#' @keywords internal
fseval.inertiaratio <-
  function (train, labels, type = c ("univariate", "multivariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (type [1] == "univariate")
    {
      return (apply (train, 2, function (v)
      {
        centers = tapply (v, labels, mean)
        center = mean (v)
        inter = (centers - center) ^2 * as.numeric (table (labels))
        inter = sum (inter)
        total = (v - center)^2
        total = sum (total)
        return (inter / total)
      }))
    }
    else
    {
      centers = apply (train, 2, function (v) tapply (v, labels, mean))
      center = matrix (apply (train, 2, mean), nrow = 1)
      inter = flexclust::dist2 (center, centers)^2 * as.numeric (table (labels))
      inter = sum (inter)
      total = flexclust::dist2 (center, train)^2
      total = sum (total)
      return (inter / total)
    }
  }

#' @keywords internal
fseval.fisher <-
  function (train, labels, type = c ("univariate", "multivariate"), ...)
  {
    if (type [1] == "univariate")
    {
      n = as.vector (table (labels))
      centers = apply (train, 2, function (v) tapply (v, labels, mean))
      center = matrix (apply (train, 2, mean), nrow = 1)
      sigmas = apply (train, 2, function (v) tapply (v, labels, stats::sd))
      return (colSums (sweep (sweep (centers, 2, center, "-")^2, 1, n, "*")) / colSums (sweep (sigmas^2, 1, n, "*")))
    }
    else
    {
      message ("Fisher score is an univariate measure")
      return (NULL)
    }
  }

#' @keywords internal
fseval.fstat <-
  function (train, labels, type = c ("univariate", "multivariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    k = nlevels (labels)
    if (type [1] == "univariate")
    {
      return (apply (train, 2, function (v)
      {
        centers = tapply (v, labels, mean)
        center = mean (v)
        inter = sum ((centers - center) ^2 * as.numeric (table (labels))) / (k - 1)
        intra = sapply (levels (labels), function (lev)
        {
          return (sum ((centers [lev] - v [labels == lev])^2))
        })
        intra = sum (intra) / (length (v) - k)
        return (inter / intra)
      }))
    }
    else
    {
      centers = apply (train, 2, function (v) tapply (v, labels, mean))
      center = matrix (apply (train, 2, mean), nrow = 1)
      inter = sum (flexclust::dist2 (center, centers)^2 * as.numeric (table (labels))) / (k - 1)
      intra = sapply (levels (labels), function (lev)
      {
        center = matrix (centers [lev, ], nrow = 1)
        return (sum (flexclust::dist2 (center, train [labels == lev, ])^2))
      })
      intra = sum (intra) / (nrow (train) - k)
      return (inter / intra)
    }
  }

#' @keywords internal
fseval.mrmr <-
  function (train, labels, type = c ("multivariate", "univariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (type [1] == "univariate")
    {
      message ("mRMR is an multivariate measure")
      return (NULL)
    }
    else
    {
      if (is.vector (train))
        train = matrix (train, ncol = 1)
      n = ncol (train)
      d = sum (apply (train, 2, function (v) mutualinformation (v, labels))) / n
      r = sum (sapply (1:n, function (i) sapply (1:n, function (j)
      {
        if (i > j)
          return (mutualinformation (train [, i], train [, j]))
        else
          return (NA)
      })), na.rm = TRUE) / (n * n)
      return (d - r)
    }
  }

#' @keywords internal
fseval.relief <-
  function (train, labels, type = c ("univariate", "multivariate"), nsamples = length (labels), k = 10, ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (type [1] == "univariate")
    {
      samples = sample (nrow (train), nsamples, replace = nsamples > nrow (train))
      dis = flexclust::dist2 (train, train)
      lev = levels (labels)
      prior = as.vector (table (labels)) / length (labels)
      names (prior) = lev
      minv = apply (train, 2, min)
      maxv = apply (train, 2, max)
      diffv = maxv - minv
      res = sapply (samples, function (s)
      {
        indices = tapply (dis [s, ], labels, function (v) order (v) [1:(k + 1)])
        res = sapply (lev, function (l)
        {
          currentclass = train [labels == l, ]
          if (is.vector (currentclass))
            currentclass = matrix (currentclass, ncol = 1)
          hnm = as.matrix (currentclass [indices [[l]], ])
          if (is.vector (hnm))
            hnm = matrix (hnm, ncol = 1)
          if (l == labels [s])
          {
            hnm = matrix (hnm [-1, ], nrow = k)
            d = abs (sweep (hnm, 2, unlist (train [s, ]), "-")) / diffv
            return (-(apply (d, 2, sum) / (nsamples * k)))
          }
          else
          {
            hnm = matrix (hnm [-1, ], nrow = k)
            d = abs (sweep (hnm, 2, unlist (train [s, ]), "-")) / diffv
            return ((prior [l] / (1 - prior [labels [s]])) * (apply (d, 2, sum) / (nsamples * k)))
          }
        })
        if (is.vector (res))
          res = matrix (res, nrow = 1)
        return (apply (res, 1, sum))
      })
      if (is.vector (res))
        res = matrix (res, nrow = 1)
      return (apply (res, 1, sum))
    }
    else
    {
      message ("Relief is an univariate measure")
      return (NULL)
    }
  }

#' @keywords internal
fseval.wrapper <-
  function (train, labels, wrapmethod, wrapeval = "accuracy", nruns = 100, ...)
  {
    return (bootstrap (methods = wrapmethod, x = train, y = labels, eval = wrapeval [1], ...))
  }

#' @keywords internal
mutualinformation <-
  function (X, Y)
  {
    uX = seq (min (X), max (X), length.out = grDevices::nclass.Sturges (X))
    coeffX = uX [2] - uX [1]
    uY = NULL
    coeffY = 1
    if (is.factor (Y))
      uY = unique (Y)
    else
    {
      uY = seq (min (Y), max (Y), length.out = grDevices::nclass.Sturges (Y))
      coeffY = uY [2] - uY [1]
    }
    n1 = length (uX)
    n2 = length (uY)
    couple = cbind.data.frame (X, Y)
    px = sapply (uX, function (x) proba (x, X))
    py = sapply (uY, function (y) proba (y, Y))
    sum (sapply (1:n1, function (i) sapply (1:n2, function (j)
    {
      xy = cbind.data.frame (uX [i], uY [j])
      pxy = proba (xy, couple)
      ppx = px [i]
      ppy = py [j]
      res = pxy * log2 (pxy / (ppx * ppy))
      if (is.nan (res))
        res = 0
      return (coeffX * coeffY * res)
    })))
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained by any classification or regression model.
#' @name predict.selection
#' @param object The classification model (of class \code{\link{cda-class}}, created by \code{\link{CDA}}).
#' @param test The test set (a \code{data.frame}).
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @param ... Other parameters.
#' @return A vector of predicted values (\code{factor}).
#' @export
#' @method predict selection
#' @seealso \code{\link{FEATURESELECTION}}, \code{\link{selection-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = FEATURESELECTION (d$train.x, d$train.y, uninb = 2, mainmethod = LDA)
#' predict (model, d$test.x)
predict.selection <-
  function (object, test, fuzzy = FALSE, ...)
  {
    test = test [, object$selection]
    if (is.vector (test))
      test = matrix (test, ncol = 1)
    return (predict (object$model, test, fuzzy, ...))
  }

#' @keywords internal
proba <-
  function (x, train)
  {
    if (is.factor (train))
      return ((as.vector (table (train)) / length (train)) [as.numeric (x)])
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    d = ncol (train)
    zs = sweep (train, 2, x, function (a, b)
    {
      if (ncol (a) == 2)
      {
        if (is.factor (a [, 2]))
        {
          select = a [, 2] == unlist (b [, 2])
          return (a [select, 1] - unlist (b [, 1]) [select])
        }
        else
          return (a - matrix (unlist (b), nrow = nrow (a)))
      }
      else
        return (a - b)
    })
    if (is.vector (zs))
      zs = matrix (zs, ncol = 1)
    h = sqrt (sum (apply (zs, 2, function (v) diff (seq (min (v), max (v), length.out = grDevices::nclass.Sturges (v)) [1:2]))^2))
    sigma = stats::cor (zs)
    invsigma = solve (sigma)
    detsigma = det (sigma)
    return (mean (apply (zs, 1, function (z)
    {
      return (exp (-(t (z) %*% invsigma %*% z) / (2 * h * h)) / ((2 * pi)^(d / 2) * h^d * sqrt (detsigma)))
    })))
  }

#' Feature selection for classification
#'
#' Select a subset of features for a classification task.
#' @name selectfeatures
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param algorithm The feature selection algorithm.
#' @param unieval The (univariate) evaluation criterion. \code{uninb}, \code{unithreshold} or \code{multieval} must be specified.
#' @param uninb The number of selected feature (univariate evaluation).
#' @param unithreshold The threshold for selecting feature (univariate evaluation).
#' @param multieval The (multivariate) evaluation criterion.
#' @param wrapmethod The classification method used for the wrapper evaluation.
#' @param keep If true, the dataset is kept in the returned result.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{FEATURESELECTION}}, \code{\link{selection-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "fstat")
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking", uninb = 2)
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
#'                 multieval = "wrapper", wrapmethod = LDA)
selectfeatures <-
  function (train,
            labels,
            algorithm = c ("ranking", "forward", "backward", "exhaustive"),
            unieval = if (algorithm [1] == "ranking") c ("fisher", "fstat", "relief", "inertiaratio") else NULL,
            uninb = NULL,
            unithreshold = NULL,
            multieval = if (algorithm [1] == "ranking") NULL else c ("mrmr", "cfs", "fstat", "inertiaratio", "wrapper"),
            wrapmethod = NULL,
            keep = FALSE,
            ...)
  {
    uni = NULL
    if (!is.null (unieval))
      uni = get (paste ("fseval.", tolower (unieval [1]), sep = ""))
    multi = NULL
    if (!is.null (multieval))
      multi = get (paste ("fseval.", tolower (multieval [1]), sep = ""))
    res = get (paste ("fs.", tolower (algorithm [1]), sep = "")) (train,
                                                                  labels,
                                                                  unieval = uni,
                                                                  uninb = uninb,
                                                                  unithreshold = unithreshold,
                                                                  multieval = multi,
                                                                  wrapmethod = wrapmethod,
                                                                  ...)
    if (!is.null (unieval))
      res$univariate = unieval [1]
    if (!is.null (multieval))
      res$multivariate = multieval [1]
    if (keep)
    {
      dataset = train [, res$selection]
      res$dataset = dataset
    }
    class (res) = "selection"
    return (res)
  }
