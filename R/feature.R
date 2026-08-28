#' Feature selection
#'
#' This class contains the result of feature selection algorithms.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{selection}}{A vector of integers indicating the selected features.}
#'   \item{\code{features}}{The names of the selected features, when the dataset had column
#'     names.}
#'   \item{\code{unieval}}{The evaluation of the features (univariate).}
#'   \item{\code{multieval}}{The evaluation of the selected features (multivariate).}
#'   \item{\code{algorithm}}{The algorithm used to select features.}
#'   \item{\code{univariate}}{The evaluation criterion (univariate).}
#'   \item{\code{nbfeatures}}{The number of features to be kept.}
#'   \item{\code{threshold}}{The threshold to decide whether a feature is kept or not.}
#'   \item{\code{multivariate}}{The evaluation criterion (multivariate).}
#'   \item{\code{dataset}}{The dataset described by the selected features only.}
#'   \item{\code{model}}{The classification model.}
#' }
#' @name selection-class
#' @seealso \code{\link{FEATURESELECTION}}, \code{\link{predict.selection}}, \code{\link{selectfeatures}}
NULL

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
#' @param mainmethod The final method used for data classification (required: either \code{mainmethod} or \code{wrapmethod} must be a valid classification/regression function, e.g. \code{LDA}, \code{NB}, ...). If a wrapper evaluation is used, the same classification method should be used.
#' @inheritParams tune.doc
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{selectfeatures}}, \code{\link{predict.selection}}, \code{\link{selection-class}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' FEATURESELECTION (iris [, -5], iris [, 5], uninb = 2, mainmethod = LDA)
#' }
FEATURESELECTION <-
  function (train,
            labels,
            algorithm = c ("ranking", "forward", "backward", "exhaustive"),
            unieval = if (algorithm [1] == "ranking") fseval.univariate () else NULL,
            uninb = NULL,
            unithreshold = NULL,
            multieval = fseval.multivariate (),
            wrapmethod = NULL,
            mainmethod = wrapmethod,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL,
            ...)
  {
    setseed (seed)
    if (!is.function (mainmethod))
      stop ("FEATURESELECTION: 'mainmethod' (or 'wrapmethod', used as its default) must be a classification/regression function, e.g. LDA, NB, KNN. Both are currently missing or invalid.")
    # tune = TRUE asks which hyperparameters a method exposes, so that performance() can
    # obtain them once and pass them back on every split. Feature selection has none of its own
    # -- the subset it keeps is not a hyperparameter, and it is chosen afresh on each split
    # anyway -- so the question goes to the method that will actually be fitted.
    if (tune)
      return (mainmethod (train, labels, tune = TRUE, ...))
    selection = selectfeatures (train, labels, algorithm, unieval, uninb, unithreshold, multieval, wrapmethod, keep = TRUE, ...)
    selection$model = mainmethod (selection$dataset, labels, methodparameters = methodparameters, ...)
    return (selection)
  }

#' @keywords internal
fs.backward <-
  function (train, labels, multieval, wrapmethod = NULL, unieval = NULL, uninb = NULL,
            unithreshold = NULL, ...)
  {
    select = rep (TRUE, ncol (train))
    besteval = multieval (train, labels, vtype = "multivariate", wrapmethod = wrapmethod, ...)
    bestselect = select
    while (sum (select) > 1)
    {
      removenext = which (select)
      eval = sapply (removenext, function (index)
      {
        nextselect = select
        nextselect [index] = FALSE
        subset = train [, nextselect]
        return (multieval (subset, labels, vtype = "multivariate", wrapmethod = wrapmethod, ...))
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
  function (train, labels, multieval, wrapmethod = NULL, unieval = NULL, uninb = NULL,
            unithreshold = NULL, ...)
  {
    p = ncol (train)
    if (p > 25)
      stop ("fs.exhaustive: ", p, " features would require evaluating 2^", p,
            " subsets, which is intractable. Please reduce the number of features, ",
            "or use algorithm = \"forward\", \"backward\", or \"ranking\" instead.")
    else if (p > 20)
      warning ("fs.exhaustive: ", p, " features require evaluating 2^", p, " = ", 2^p,
               " subsets; this may take a long time. Consider algorithm = \"forward\" ",
               "or \"backward\" for larger feature sets.")
    indices = 1:(2^ncol (train) - 1)
    eval = sapply (indices, function (index)
    {
      select = as.logical (as.integer (intToBits (index) [1:ncol (train)]))
      subset = train [, select]
      return (multieval (subset, labels, vtype = "multivariate", wrapmethod = wrapmethod, ...))
    })
    res = as.logical (as.integer (intToBits (which.max (eval)) [1:ncol (train)]))
    res = list (selection = which (res), multieval = max (eval), algorithm = "exhaustive")
    return (res)
  }

#' @keywords internal
fs.forward <-
  function (train, labels, multieval, wrapmethod = NULL, unieval = NULL, uninb = NULL,
            unithreshold = NULL, ...)
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
        return (multieval (subset, labels, vtype = "multivariate", wrapmethod = wrapmethod, ...))
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
    eval = unieval (train, labels, vtype = "univariate")
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
    else if (!is.null (multieval))
    {
      size = 1:ncol (train)
      features = order (eval, decreasing = TRUE)
      meval = sapply (size, function (index)
      {
        subset = train [, features [1:index]]
        return (multieval (subset, labels, vtype = "multivariate", wrapmethod = wrapmethod, ...))
      })
      selection = sort (features [1:which.max (meval)])
      res = list (selection = selection, unieval = eval, multieval = max (meval), algorithm = "ranking")
    }
    else
      stop ("selectfeatures: with algorithm = \"ranking\", one of 'uninb' (how many features ",
            "to keep), 'unithreshold' (the score above which a feature is kept) or ",
            "'multieval' (a multivariate criterion used to choose how many of the ranked ",
            "features to keep) must be given. None of them was.")
    return (res)
  }

#' @keywords internal
fseval.cfs <-
  function (train, labels, vtype = c ("multivariate", "univariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (vtype [1] == "univariate")
    {
      message ("CFS is a multivariate measure")
      return (NULL)
    }
    else
    {
      # Hall's merit: k * mean (r_cf) / sqrt (k + k (k - 1) mean (r_ff)), where r_cf is the
      # correlation between a feature and the class and r_ff that between two features.
      #
      # r_cf is the correlation ratio -- the square root of the between-class share of the
      # variance, which fseval.inertiaratio() already computes. It is the usual generalisation
      # of |r| to a numeric feature and a categorical class: bounded by [0, 1], exactly |r| on
      # two classes, and deterministic.
      k = ncol (train)
      rcf = mean (sqrt (fseval.inertiaratio (train, labels, vtype = "univariate")))
      rff = 0
      if (k > 1)
      {
        # A constant variable has no standard deviation, so cor() returns NA for every pair it
        # belongs to and the whole merit came back NA. It is not redundant with anything
        # either, which is what a correlation of 0 says.
        r = suppressWarnings (stats::cor (train)) [lower.tri (diag (k))]
        r [is.na (r)] = 0
        rff = mean (abs (r))
      }
      return (k * rcf / sqrt (k + (k * (k - 1) * rff)))
    }
  }

#' @keywords internal
fseval.inertiaratio <-
  function (train, labels, vtype = c ("univariate", "multivariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (vtype [1] == "univariate")
    {
      return (apply (train, 2, function (v)
      {
        centers = tapply (v, labels, mean)
        center = mean (v)
        inter = (centers - center) ^2 * as.numeric (table (labels))
        inter = sum (inter)
        total = (v - center)^2
        total = sum (total)
        # A constant variable has no variance to share out; it separates nothing, so its ratio
        # is 0 rather than 0 / 0.
        if (total == 0)
          return (0)
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
      if (total == 0)
        return (0)
      return (inter / total)
    }
  }

#' @keywords internal
fseval.fisher <-
  function (train, labels, vtype = c ("univariate", "multivariate"), ...)
  {
    if (vtype [1] == "univariate")
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
  function (train, labels, vtype = c ("univariate", "multivariate"), ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    k = nlevels (labels)
    if (vtype [1] == "univariate")
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
  function (train, labels, vtype = c ("multivariate", "univariate"), micache = NULL, k = 3, ...)
  {
    if (vtype [1] == "univariate")
    {
      message ("mRMR is a multivariate measure")
      return (NULL)
    }
    # A single column arrives as a bare vector. It is wrapped in a one-column data.frame
    # rather than a matrix so that a factor stays a factor: apply() and as.matrix() turn a
    # mixed data.frame into a character matrix, which would have every numeric variable
    # estimated as if it were categorical.
    if (is.null (dim (train)))
      train = data.frame (X = train)
    p = ncol (train)
    column = if (is.data.frame (train)) function (j) train [[j]] else function (j) train [, j]
    # Feature selection evaluates the same pair of variables again and again -- backward
    # elimination on q features asks for O(q^4) pair estimates, of which only q(q-1)/2 are
    # distinct. selectfeatures() therefore hands this function an environment in which each
    # estimate is kept, keyed by the names of the variables involved. Without such a cache
    # (direct call, or a matrix whose columns are not uniquely named) every estimate is simply
    # recomputed, as before.
    names.ok = (!is.null (colnames (train))) && (!anyNA (colnames (train))) &&
               all (nzchar (colnames (train))) && (anyDuplicated (colnames (train)) == 0)
    cached = (!is.null (micache)) && names.ok
    nm = colnames (train)
    mi = function (key, a, b)
    {
      if (!cached)
        return (mutualinformation (a, b, k))
      key = paste (key, k, sep = "\r")
      hit = micache [[key]]
      if (!is.null (hit))
        return (hit)
      value = mutualinformation (a, b, k)
      assign (key, value, envir = micache)
      return (value)
    }
    relevance = sum (sapply (1:p, function (j)
      mi (paste ("\r\ry", nm [j], sep = "\r"), column (j), labels))) / p
    redundancy = 0
    if (p > 1)
    {
      pairs = utils::combn (p, 2)
      redundancy = sum (apply (pairs, 2, function (idx)
        mi (paste (sort (nm [idx]), collapse = "\r"), column (idx [1]), column (idx [2])))) /
        (p * p)
    }
    return (relevance - redundancy)
  }

#' @keywords internal
# The importance a random forest gives each variable: how much accuracy the forest loses when
# that variable's values are shuffled, averaged over the trees. Unlike the other univariate
# criteria it does not look at a variable on its own -- the forest weighs each one in the
# presence of the others -- but it answers the same question and ranks in the same direction,
# so it plugs in as one of them.
fseval.randomforest <-
  function (train, labels, vtype = c ("univariate", "multivariate"), ntree = 500, seed = NULL, ...)
  {
    if (vtype [1] != "univariate")
    {
      message ("Random forest importance is an univariate measure")
      return (NULL)
    }
    if (!requireNamespace ("randomForest", quietly = TRUE))
      stop ("selectfeatures: unieval = \"randomforest\" needs the 'randomForest' package. ",
            "Please install it, or use another criterion.")
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    train = as.data.frame (train)
    setseed (seed)
    forest = randomForest::randomForest (x = train, y = labels, ntree = ntree,
                                         importance = TRUE)
    # type = 1 is the permutation importance -- the mean decrease in accuracy, or in MSE for a
    # numeric target. type = 2, the node impurity, favours the variables that take many
    # distinct values whether or not they are informative.
    res = randomForest::importance (forest, type = 1) [, 1]
    names (res) = colnames (train)
    return (res)
  }

#' @keywords internal
fseval.relief <-
  function (train, labels, vtype = c ("univariate", "multivariate"),
            nsamples = length (labels), k = 10, ...)
  {
    if (is.vector (train))
      train = matrix (train, ncol = 1)
    if (vtype [1] != "univariate")
    {
      message ("Relief is an univariate measure")
      return (NULL)
    }
    train = as.matrix (train)
    labels = factor (labels)
    n = nrow (train)
    p = ncol (train)
    lev = levels (labels)
    prior = as.vector (table (labels)) / n
    names (prior) = lev
    # A constant variable would divide by zero; it separates nothing, so its weight is 0
    # whatever the denominator.
    diffv = apply (train, 2, function (v) diff (range (v)))
    diffv [diffv == 0] = 1
    samples = sample (n, nsamples, replace = nsamples > n)
    dis = flexclust::dist2 (train, train)
    w = numeric (p)
    names (w) = colnames (train)
    for (l in lev)
    {
      cls = which (labels == l)
      same = labels [samples] == l
      # The sampled instances are treated in two groups -- those of class 'l' and the others
      # -- because they do not have the same neighbours available: an instance of class 'l'
      # is its own nearest neighbour there and that one is dropped. Within a group everything
      # is done in one matrix operation. 'kl' keeps the number of neighbours within reach of
      # the smallest class, and the normalisation divides each feature by its own range.
      for (issame in c (TRUE, FALSE))
      {
        take = which (same == issame)
        kl = min (k, length (cls) - as.integer (issame))
        if ((length (take) == 0) || (kl < 1))
          next
        nn = knn.indices (dis [samples [take], cls, drop = FALSE], kl, drop.first = issame)
        # Rows of 'train' to compare each sampled instance with, and the coefficient the
        # comparison enters the weight with: negative for a hit, and for a miss the prior of
        # the class it comes from, renormalised over the classes that are not the instance's.
        rows = cls [as.vector (nn)]
        coefficient = if (issame) rep (-1, length (take))
                      else prior [l] / (1 - prior [as.character (labels [samples [take]])])
        d = abs (train [rows, , drop = FALSE] -
                 train [rep (samples [take], times = kl), , drop = FALSE])
        w = w + colSums (d * rep (coefficient, times = kl)) / diffv / (nsamples * kl)
      }
    }
    return (w)
  }

#' @keywords internal
# The indices of the 'k' smallest values of every row of a distance matrix, as a
# (nrow x k) matrix. 'drop.first' says, per row, whether the very nearest one is to be
# skipped -- which it is when the row's own observation is among the candidates.
#
# order() sorts the whole row, which costs O(n log n) per row and dominated relief() on
# anything but small data. Only the first few values are needed, so the row is cut at its
# (k + 1)-th smallest value first: everything strictly below it is kept outright, and the ties
# at that value make up the difference.
knn.indices <-
  function (d, k, drop.first)
  {
    n = nrow (d)
    drop.first = rep (drop.first, length.out = n)
    # One column of the answer per pass: max.col() finds the nearest remaining neighbour of
    # every row at once, and that neighbour is then masked out. k + 1 passes over the whole
    # matrix, entirely inside R's C code, rather than one sort per row.
    m = min (k + as.integer (any (drop.first)), ncol (d))
    neg = -as.matrix (d)
    rows = seq_len (n)
    res = matrix (0L, nrow = n, ncol = m)
    for (j in seq_len (m))
    {
      idx = max.col (neg, ties.method = "first")
      res [, j] = idx
      neg [cbind (rows, idx)] = -Inf
    }
    out = matrix (0L, nrow = n, ncol = k)
    plain = !drop.first
    if (any (plain))
      out [plain, ] = res [plain, seq_len (k), drop = FALSE]
    if (any (drop.first))
      out [drop.first, ] = res [drop.first, 1L + seq_len (k), drop = FALSE]
    return (out)
  }

#' @keywords internal
fseval.wrapper <-
  function (train, labels, wrapmethod, wrapeval = "accuracy", nruns = 10, ...)
  {
    if (is.null (wrapmethod))
      stop ("selectfeatures: multieval = \"wrapper\" evaluates a feature subset by actually ",
            "fitting a model on it, so 'wrapmethod' must name the classification or ",
            "regression function to use (e.g. wrapmethod = LDA). It is missing.")
    # 'nruns' was declared and then never passed on, so it said 100 and performance() did 10.
    res = performance (methods = wrapmethod, train.x = train, train.y = labels,
                       type = "evaluation", protocol = "bootstrap", eval = wrapeval [1],
                       nruns = nruns, ...)
    return (res)
  }

#' @keywords internal
# Should this variable be treated as discrete? Factors, characters and logicals always are;
# a numeric variable is treated as discrete when it takes very few distinct values, since the
# nearest-neighbour estimators below assume a continuous distribution (no ties).
mi.isdiscrete <-
  function (v, maxlevels = 10)
    is.factor (v) || is.character (v) || is.logical (v) || (length (unique (v)) <= maxlevels)

#' @keywords internal
# Breaks the ties of a continuous variable with noise several orders of magnitude below its
# own scale, as Kraskov et al. recommend: the k-nearest-neighbour estimators are undefined
# when many points share a coordinate. The draw is deterministic (fixed seed) and the user's
# random number generator is left exactly as it was found.
mi.untie <-
  function (v)
  {
    if (!anyDuplicated (v))
      return (v)
    s = stats::sd (v)
    if ((!is.finite (s)) || (s == 0))
      return (v)
    hasseed = exists (".Random.seed", envir = globalenv (), inherits = FALSE)
    if (hasseed)
    {
      saved = get (".Random.seed", envir = globalenv (), inherits = FALSE)
      on.exit (assign (".Random.seed", saved, envir = globalenv ()))
    }
    else
      on.exit (suppressWarnings (rm (".Random.seed", envir = globalenv ())))
    set.seed (1)
    return (v + stats::runif (length (v), -1e-10 * s, 1e-10 * s))
  }

#' @keywords internal
# Brings a variable to unit variance, leaving a constant one alone.
mi.scale <-
  function (v)
  {
    s = stats::sd (v)
    if ((!is.finite (s)) || (s == 0))
      return (v)
    return (v / s)
  }

#' @keywords internal
# The k-th smallest value of every row of a matrix (the point's own zero distance included,
# hence the k + 1).
mi.kth <-
  function (m, k) apply (m, 1, function (r) sort.int (r, partial = k + 1) [k + 1])

#' @keywords internal
# Mutual information between two continuous variables: the k-nearest-neighbour estimator of
# Kraskov, Stogbauer & Grassberger (2004), estimator (1). Distances are taken in the maximum
# norm, which is what makes the volume terms of the two marginal spaces cancel.
# The n x n distance matrices are built one block of rows at a time so that memory stays
# bounded on large samples.
mi.ksg <-
  function (x, y, k = 3, chunk = 2048)
  {
    n = length (x)
    k = min (k, n - 1)
    if (k < 1)
      return (0)
    # Both variables are brought to unit variance first. A mutual information does not depend
    # on the units the variables are measured in, but the maximum norm below mixes the two
    # coordinates, so the estimate does unless they are on a common scale: on a bivariate
    # normal whose true value is 0.322 bit, multiplying one variable by 1000 turned the
    # estimate into 0.005 without this. Scaling is free -- an affine change of variable leaves
    # the mutual information unchanged.
    x = mi.scale (x)
    y = mi.scale (y)
    nx = integer (n)
    ny = integer (n)
    for (s in seq (1, n, by = chunk))
    {
      e = min (s + chunk - 1, n)
      dx = abs (outer (x [s:e], x, "-"))
      dy = abs (outer (y [s:e], y, "-"))
      eps = mi.kth (pmax (dx, dy), k)
      nx [s:e] = rowSums (dx < eps) - 1L
      ny [s:e] = rowSums (dy < eps) - 1L
    }
    return ((digamma (k) + digamma (n) - mean (digamma (nx + 1) + digamma (ny + 1))) / log (2))
  }

#' @keywords internal
# Mutual information between a continuous and a discrete variable: the k-nearest-neighbour
# estimator of Ross (2014). For each observation, the k-th nearest neighbour is looked for
# among the observations of its own class only, and the resulting radius is then used to
# count neighbours in the whole sample. 'k' is reduced for the classes that are too small to
# provide k neighbours.
mi.ross <-
  function (x, y, k = 3, chunk = 2048)
  {
    y = factor (y)
    n = length (x)
    d = numeric (n)
    kk = integer (n)
    for (l in levels (y))
    {
      idx = which (y == l)
      kl = min (k, length (idx) - 1)
      if (kl < 1)
        next
      d [idx] = mi.kth (abs (outer (x [idx], x [idx], "-")), kl)
      kk [idx] = kl
    }
    keep = kk > 0
    if (!any (keep))
      return (0)
    m = integer (n)
    for (s in seq (1, n, by = chunk))
    {
      e = min (s + chunk - 1, n)
      m [s:e] = rowSums (abs (outer (x [s:e], x, "-")) <= d [s:e]) - 1L
    }
    ny = as.vector (table (y)) [as.integer (y)]
    return ((digamma (n) - mean (digamma (ny [keep])) +
             mean (digamma (kk [keep])) - mean (digamma (pmax (m [keep], 1)))) / log (2))
  }

#' @keywords internal
# Mutual information between two discrete variables: the plug-in estimator on the observed
# joint distribution, corrected for its (always positive) small-sample bias by the
# Miller-Madow correction. Truncated at zero, since a mutual information cannot be negative.
mi.discrete <-
  function (x, y)
  {
    tt = table (factor (x), factor (y))
    n = sum (tt)
    if (n == 0)
      return (0)
    p = tt / n
    joint = outer (rowSums (p), colSums (p))
    keep = p > 0
    raw = sum (p [keep] * log2 (p [keep] / joint [keep]))
    correction = (sum (tt > 0) - sum (rowSums (tt) > 0) - sum (colSums (tt) > 0) + 1) /
                 (2 * n * log (2))
    return (max (0, raw - correction))
  }

#' @keywords internal
# Mutual information between two variables, in bits. This is the quantity mRMR is defined on
# (see fseval.mrmr), and the estimator is chosen according to the nature of the two variables:
#
#   continuous / continuous  Kraskov, Stogbauer & Grassberger (2004), estimator (1)
#   continuous / discrete    Ross (2014)
#   discrete / discrete      plug-in with the Miller-Madow bias correction
#
# These are the reference estimators for the problem, and are an order of magnitude more
# accurate than a binned one: on a bivariate normal with 150 observations and a true mutual
# information of 0.322 bit, a histogram returns 0.536 and this returns 0.341. The cost is
# O(n^2) per pair of variables rather than O(n); see fseval.mrmr, which caches its estimates so
# that each pair is only ever computed once.
mutualinformation <-
  function (X, Y, k = 3)
  {
    if (length (X) != length (Y))
      stop ("mutualinformation: the two variables must have the same length (", length (X),
            " and ", length (Y), ").")
    dx = mi.isdiscrete (X)
    dy = mi.isdiscrete (Y)
    if (dx && dy)
      return (mi.discrete (X, Y))
    if (dx)
      return (max (0, mi.ross (mi.untie (as.numeric (Y)), X, k)))
    if (dy)
      return (max (0, mi.ross (mi.untie (as.numeric (X)), Y, k)))
    return (max (0, mi.ksg (mi.untie (as.numeric (X)), mi.untie (as.numeric (Y)), k)))
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
#' \dontrun{
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = FEATURESELECTION (d$train.x, d$train.y, uninb = 2, mainmethod = LDA)
#' predict (model, d$test.x)
#' }
predict.selection <-
  function (object, test, fuzzy = FALSE, ...)
  {
    # By name whenever the selection recorded them: indexing by position picks the wrong
    # variables, silently, as soon as 'test' does not carry its columns in the order the
    # training set had.
    columns = object$selection
    if ((!is.null (object$features)) && (!is.null (colnames (test))))
    {
      unknown = setdiff (object$features, colnames (test))
      if (length (unknown) > 0)
        stop ("predict.selection: 'test' does not have the selected variable(s): ",
              paste (unknown, collapse = ", "), ".")
      columns = object$features
    }
    test = test [, columns, drop = FALSE]
    res = predict (object$model, test, fuzzy, ...)
    return (res)
  }

#' Plot a feature selection
#'
#' Draws the score every variable obtained, the ones that were kept apart from the ones that
#' were not. Only a selection made by ranking the variables (\code{algorithm = "ranking"}) can
#' be drawn this way: it is the only one that scores them one by one, where the other three
#' score whole subsets.
#' @name plot.selection
#' @param x The selection (object of class \code{\link{selection-class}}, created by
#' \code{\link{selectfeatures}}).
#' @param horiz Whether the bars are drawn horizontally, which leaves room for long variable
#' names.
#' @param legendpos Position of the legend.
#' @param ... Other parameters, passed to \code{\link[graphics]{barplot}}.
#' @export
#' @method plot selection
#' @seealso \code{\link{selectfeatures}}, \code{\link{selection-class}},
#' \code{\link{print.selection}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' # How useful a random forest finds each variable, and the two it would keep
#' selection = selectfeatures (iris [, -5], iris [, 5], unieval = "randomforest", uninb = 2)
#' selection
#' plot (selection)
#' }
plot.selection <-
  function (x, horiz = TRUE, legendpos = "bottomright", ...)
  {
    scores = x$unieval
    if (is.null (scores) || (length (scores) == 0))
      stop ("plot.selection: this selection carries no score per variable. Only ",
            "algorithm = \"ranking\" scores the variables one by one; \"forward\", ",
            "\"backward\" and \"exhaustive\" score whole subsets, which cannot be drawn ",
            "this way.")
    if (is.null (names (scores)))
      names (scores) = paste ("V", seq_along (scores))
    # Horizontally, barplot() draws the first bar at the bottom, so the best score has to come
    # last for the picture to read from the top down.
    order = order (scores, decreasing = !horiz)
    kept = order %in% x$selection
    scores = scores [order]
    colours = ifelse (kept, 2, "grey75")
    if (horiz)
    {
      extra = max (graphics::strwidth (names (scores), units = "figure") * 30)
      opar = graphics::par (mar = graphics::par ("mar") + c (0, extra, 0, 0))
      on.exit (graphics::par (opar))
    }
    graphics::barplot (scores, horiz = horiz, col = colours, border = NA, las = 1,
                       xlab = if (horiz) x$univariate else "",
                       ylab = if (horiz) "" else x$univariate, ...)
    graphics::legend (legendpos, c ("kept", "dropped"), fill = c (2, "grey75"), bty = "n")
  }

#' @keywords internal
# Lookup table for the fs.* dispatch used by selectfeatures(). Built lazily (as a function
# rather than a static list) so it does not depend on file collation order.
#
# Every fs.* function takes the *same* argument list, whether it uses all of it or not, because
# selectfeatures() passes all of it by name -- an argument left to fall into '...' would travel
# on to the evaluation criterion and, through fseval.wrapper(), into the learning method.
fs.functions <-
  function ()
  {
    list (backward = fs.backward,
          exhaustive = fs.exhaustive,
          forward = fs.forward,
          ranking = fs.ranking)
  }

#' @keywords internal
# The criteria that make sense on a single variable, and those that need a subset. Kept next
# to the dispatch table so that adding a criterion means touching one place.
fseval.univariate <-
  function () c ("fisher", "fstat", "relief", "inertiaratio", "randomforest")

#' @keywords internal
fseval.multivariate <-
  function () c ("mrmr", "cfs", "fstat", "inertiaratio", "wrapper")

#' @keywords internal
# Lookup table for the fseval.* dispatch used by selectfeatures(). Built lazily (as a
# function rather than a static list) so it does not depend on file collation order.
fseval.functions <-
  function ()
  {
    list (cfs = fseval.cfs,
          fisher = fseval.fisher,
          fstat = fseval.fstat,
          inertiaratio = fseval.inertiaratio,
          mrmr = fseval.mrmr,
          randomforest = fseval.randomforest,
          relief = fseval.relief,
          wrapper = fseval.wrapper)
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
#' \dontrun{
#' require (datasets)
#' data (iris)
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "forward", multieval = "fstat")
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking", uninb = 2)
#' selectfeatures (iris [, -5], iris [, 5], algorithm = "ranking",
#'                 multieval = "wrapper", wrapmethod = LDA)
#' }
selectfeatures <-
  function (train,
            labels,
            algorithm = c ("ranking", "forward", "backward", "exhaustive"),
            unieval = if (algorithm [1] == "ranking") fseval.univariate () else NULL,
            uninb = NULL,
            unithreshold = NULL,
            multieval = fseval.multivariate (),
            wrapmethod = NULL,
            keep = FALSE,
            ...)
  {
    fs.funs = fs.functions ()
    fseval.funs = fseval.functions ()
    algorithm = match.arg (algorithm [1], names (fs.funs))
    # Derived from the dispatch table, so that a criterion is declared in a single place.
    uni = NULL
    if (!is.null (unieval))
    {
      unieval = match.arg (unieval [1], fseval.univariate ())
      uni = fseval.funs [[tolower (unieval)]]
    }
    multi = NULL
    if (!is.null (multieval))
    {
      multieval = match.arg (multieval [1], fseval.multivariate ())
      multi = fseval.funs [[tolower (multieval)]]
      if (tolower (multieval) == "mrmr")
      {
        # One cache per call to selectfeatures(), handed to fseval.mrmr() through a closure
        # rather than through '...': the other criteria forward '...' to a learning method
        # (fseval.wrapper() does), and an unexpected 'micache' argument would land there.
        micache = new.env (parent = emptyenv ())
        mrmr = multi
        multi = function (train, labels, ...) mrmr (train, labels, ..., micache = micache)
      }
    }
    res = fs.funs [[algorithm]] (train,
                                 labels,
                                 unieval = uni,
                                 uninb = uninb,
                                 unithreshold = unithreshold,
                                 multieval = multi,
                                 wrapmethod = wrapmethod,
                                 ...)
    if (is.null (res) || is.null (res$selection) || (length (res$selection) == 0))
      stop ("selectfeatures: algorithm \"", algorithm, "\" did not select any feature. ",
            "Check the criteria and thresholds passed to it.")
    if (!is.null (unieval))
      res$univariate = unieval
    # Only when it was actually used: ranking by a univariate criterion with 'uninb' or
    # 'unithreshold' never calls the multivariate one, and print.selection() announced a
    # criterion that had played no part.
    if ((!is.null (multieval)) && ((algorithm != "ranking") || (!is.null (res$multieval))))
      res$multivariate = multieval
    # The names of the retained variables, not only their positions: print.selection() shows
    # them, and they are what one wants to read back out of the result.
    if (!is.null (colnames (train)))
      res$features = colnames (train) [res$selection]
    if (keep)
      res$dataset = train [, res$selection, drop = FALSE]
    class (res) = "selection"
    return (res)
  }
