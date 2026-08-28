fdm2id.globals <- new.env (emptyenv ())
fdm2id.globals$export <- TRUE

#' Shared documentation of the arguments every learning method takes
#'
#' This function is never called: it holds the canonical documentation of the four arguments
#' that close the signature of every classification and regression method of the package,
#' shared through \code{@@inheritParams} rather than repeated in some thirty places. A method
#' whose own \code{@@param} says something more specific keeps it.
#'
#' Every learning method of the package ends on the same four arguments, in the same order:
#' \code{tune}, \code{methodparameters}, \code{graph}, \code{seed}. That is what lets
#' \code{\link{performance}} take any of them without knowing which, and what lets one method
#' be replaced by another in a script without rewriting the call.
#' @name tune.doc
#' @param tune If true, the function returns parameters instead of a classification model.
#' @param nfolds The number of folds of the cross-validation a method runs to choose its
#' hyperparameters. Only used when there is something to choose, i.e. when one of them is given
#' as a vector. Lower it to fit faster, at the cost of a noisier choice.
#' @param methodparameters Pre-tuned parameters, as returned by the same method called with
#' \code{tune = TRUE}. \code{\link{performance}} obtains them once and passes them back when
#' fitting, so that the tuning is not redone on every split. A method with nothing to tune
#' returns an empty object, which leaves its defaults untouched.
#' @param graph Whether the method draws the graphic that goes with its tuning (the
#' cross-validation curve, typically). Methods that have no such graphic accept the argument
#' and ignore it.
#' @param seed A specified seed for random number generation, so that two runs on the same
#' data give the same model. Every learning method accepts it, so that it can be set the same
#' way whatever the method; the deterministic ones simply have nothing to draw and give the
#' same model with or without it.
#' @keywords internal
# The formals must list every argument documented above, or R CMD check reports the extra ones
# as "documented arguments not in \usage".
tune.doc <-
  function (tune, methodparameters, graph, seed, nfolds) NULL

#' @keywords internal
# The resampling scheme every method searching a grid of hyperparameters hands to e1071's
# tune.*() functions -- a plain nfolds-fold cross-validation, which is both e1071's own default
# and the usual choice for model selection.
tune.scheme <-
  function (nfolds = 10)
  {
    if ((!is.numeric (nfolds)) || (length (nfolds) != 1) || is.na (nfolds) || (nfolds < 2))
      stop ("'nfolds' must be a single whole number >= 2 (the number of folds of the ",
            "cross-validation used to choose the hyperparameters). Got: ", nfolds, ".")
    return (e1071::tune.control (sampling = "cross", cross = nfolds))
  }

#' @keywords internal
# Cluster labels arrive either as integer codes (0 marking the observations a density method
# left as noise) or as a factor / character vector -- a ground truth, the output of a
# classifier. Everything downstream expects the codes (min(), 1 + clusters, colour indices), so
# categorical input is converted here and its labels kept for the legend.
cluster.codes <-
  function (clusters)
  {
    names = NULL
    if (is.factor (clusters) || is.character (clusters) || is.logical (clusters))
    {
      clusters = factor (clusters)
      names = levels (clusters)
      clusters = as.numeric (clusters)
    }
    return (list (codes = clusters, names = names))
  }

#' @keywords internal
# The legend of a clustering plot, and the colour of each entry: the names of the classes when
# the input was categorical, "Cluster i" otherwise, and "Noise" for cluster 0. Colours follow
# the 1 + code convention used everywhere else, so noise is black.
cluster.legend <-
  function (codes, names = NULL)
  {
    values = sort (unique (codes))
    labels = paste ("Cluster", values)
    if (!is.null (names))
      labels = names [values]
    else if (min (codes) == 0)
      labels = c ("Noise", paste ("Cluster", values [-1]))
    return (list (labels = labels, col = 1 + values))
  }

#' @keywords internal
# Seeds the random number generator, but *only* when a seed was actually given.
#
# set.seed (NULL) is not a no-op: it re-initialises the generator from the current time and the
# process id, silently throwing away a seed the caller had just set, so that
#
#   set.seed (42); splitdata (iris, 5)
#
# would not be reproducible. This helper leaves the generator untouched when 'seed' is NULL --
# or empty, which happens with expressions such as 'seed + i' -- so an ambient seed is honoured.
setseed <-
  function (seed)
  {
    if ((!is.null (seed)) && (length (seed) > 0) && (!is.na (seed [1])))
      set.seed (seed [1])
    invisible (NULL)
  }

#' @keywords internal
# Stops with an explicit message when a set of predictors handed to a method that only accepts
# numbers -- the support vector machines -- carries a categorical column.
#
# e1071::svm() does not check: it coerces the factor into the design matrix and fails much
# later, inside its own scaling, with "missing value where TRUE/FALSE needed". That message
# names neither the function nor the column, and the practical sessions of the course use this
# very failure to make the point that an SVM needs numbers, so it has to be readable.
check.numeric.predictors <-
  function (d, fname)
  {
    if (is.vector (d) || is.factor (d))
      d = data.frame (X = d)
    bad = which (!sapply (as.data.frame (d), is.numeric))
    if (length (bad) > 0)
      stop (fname, ": the predictors must all be numeric, but ",
            ifelse (length (bad) == 1, "column ", "columns "),
            paste (names (bad), collapse = ", "),
            ifelse (length (bad) == 1, " is ", " are "),
            "categorical. A support vector machine works on distances between observations, ",
            "which a factor does not define. Either keep the numeric columns only, or recode ",
            "each factor as indicator variables -- one numeric column per modality.")
    invisible (NULL)
  }

#' @keywords internal
# Resolves a 'target' argument given either as a column index or as a column name, so that
# every function taking one accepts both.
column.index <-
  function (dataset, target, fname)
  {
    if (!is.character (target))
      return (target)
    index = which (colnames (dataset) == target)
    if (length (index) == 0)
      stop (fname, ": 'target' (\"", target, "\") does not match any column name of 'dataset'. ",
            "Available columns: ", paste (colnames (dataset), collapse = ", "))
    return (index)
  }

#' @keywords internal
# Adds an alpha channel to a vector of colours (names, hex strings or palette indices).
# col2rgb() is vectorised, which matters: plotdata() calls this with one colour per observation.
addalpha <-
  function (colors, a = 32)
  {
    rgb = grDevices::col2rgb (colors)
    return (grDevices::rgb (rgb [1, ], rgb [2, ], rgb [3, ], a, maxColorValue = 255))
  }

#' Duplicate and add noise to a dataset
#'
#' This function is a data augmentation technique. It duplicates rows and add gaussian noise to the duplicates.
#' @name augmentation
#' @param dataset The dataset to be split (\code{data.frame} or \code{matrix}).
#' @param target The column index (numeric) or column name (character) of the target variable (class label or response variable).
#' @param n The scaling factor (as an integer value): the output contains \code{n} times the original dataset (the original rows, plus \code{n - 1} noisy copies).
#' @param sigma The baseline variance for the noise generation.
#' @param seed A specified seed for random number generation.
#' @return An augmented dataset.
#' @export
#' @examples
#' require (datasets)
#' data (iris)
#' d = augmentation (iris, 5)
#' summary (iris)
#' summary (d)
#' # 'target' can also be given as a column name
#' d = augmentation (iris, "Species")
augmentation <-
  function (dataset, target, n = 5, sigma = .1, seed = NULL)
  {
    target = column.index (dataset, target, "augmentation")
    if ((!is.numeric (n)) || (length (n) != 1) || is.na (n) || (n < 1) || (n != round (n)))
      stop ("augmentation: 'n' must be a single whole number >= 1 (it is the scaling factor: the ",
            "output contains 'n' times the original dataset). Got: ", n)
    # Gaussian noise only means something on a quantitative variable. apply() over a mixed
    # data.frame coerces everything to character first, so sd() returned NA and every value of
    # every copy came out NA. The qualitative variables are copied as they are.
    predictors = setdiff (seq_len (ncol (dataset)), target)
    numeric = predictors [sapply (dataset [, predictors, drop = FALSE], is.numeric)]
    if (length (numeric) == 0)
      stop ("augmentation: none of the predictors is numeric, so there is nothing to add ",
            "gaussian noise to.")
    if (length (numeric) < length (predictors))
      message ("augmentation: ", length (predictors) - length (numeric), " qualitative ",
               "variable(s) copied unchanged -- noise applies to the numeric ones only.")
    std = sapply (dataset [, numeric, drop = FALSE], stats::sd) * sigma
    copy = dataset [rep (seq_len (nrow (dataset)), n - 1), ]
    setseed (seed)
    noise = mapply (function (mu, s) stats::rnorm (nrow (copy), mu, s), mu = 0, s = std)
    copy [, numeric] = copy [, numeric] + noise
    res = rbind (dataset, copy)
    return (res)
  }

#' Correlated variables
#'
#' Return the list of correlated variables
#' @name correlated
#' @param d A data matrix.
#' @param threshold The threshold on the (absolute) Pearson coefficient. If NULL, return the most correlated variables.
#' @return The list of correlated variables (as a matrix of column names).
#' @seealso \code{\link[stats]{cor}}
#' @export
#' @examples
#' data (iris)
#' correlated (iris)
correlated <-
  function (d, threshold = 0.8)
  {
    factors = NULL
    if (is.factor (d))
      factors = TRUE
    else if (is.vector (d))
      factors = FALSE
    else
      factors = sapply (as.data.frame (d), is.factor)
    if (sum (factors) > 0)
      d = d [, !factors, drop = FALSE]
    if (is.null (ncol (d)) || (ncol (d) < 2))
      stop ("correlated: at least two numeric variables are needed to look for correlated ",
            "pairs; 'd' has ", if (is.null (ncol (d))) 1 else ncol (d),
            " once the qualitative variables have been removed.")
    cm = stats::cor (d)
    n = colnames (d)
    l = length (n)
    res = NULL
    if (is.null (threshold))
    {
      # The largest |r| off the diagonal, computed on a copy so that 'cm' keeps the *signed*
      # coefficients the result reports.
      offdiag = abs (cm)
      diag (offdiag) = 0
      threshold = max (offdiag)
    }
    idx = which (lower.tri (cm) & (abs (cm) >= threshold), arr.ind = TRUE)
    if (nrow (idx) == 0)
    {
      # No pair reaches the threshold: an empty table of the right shape rather than an error,
      # since raising the threshold is the first thing one does when exploring a correlation
      # structure.
      res = data.frame (character (0), character (0), numeric (0), stringsAsFactors = FALSE)
      colnames (res) = c ("Var. 1", "Var. 2", "r")
      return (res)
    }
    val = cm [idx]
    # Order each (row, column) pair so the smaller index comes first, then sort the pairs --
    # carrying 'val' along with them, so that every coefficient stays attached to its pair.
    idx = t (apply (idx, 1, sort))
    o = order (idx [, 1], idx [, 2])
    idx = idx [o, , drop = FALSE]
    val = val [o]
    res = cbind.data.frame (matrix (n [idx], ncol = 2), r = val, stringsAsFactors = FALSE)
    colnames (res) = c ("Var. 1", "Var. 2", "r")
    # By decreasing *strength* of the correlation: ordering on the signed coefficient sent
    # the strongest negative pairs to the bottom of a table meant to show the most correlated
    # variables first.
    res = res [order (-abs (val)), , drop = FALSE]
    rownames (res) = seq_len (nrow (res))
    return (res)
  }

#' Close a graphics device
#'
#' Close the graphics device driver
#' @name closegraphics
#' @param export If given, explicitly overrides the global export toggle set by \code{\link{toggleexport}} for
#' this call only (\code{TRUE} closes the device, \code{FALSE} is a no-op). By default, the global toggle is used,
#' so existing code is unaffected.
#' @seealso \code{\link{exportgraphics}}, \code{\link{toggleexport}}, \code{\link[grDevices]{dev.off}}
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' # Explicit override, ignoring the global toggle:
#' closegraphics (export = TRUE)
#' }
closegraphics <-
  function (export = fdm2id.globals$export)
  {
    if (export)
      grDevices::dev.off ()
  }

#' @keywords internal
# Maps file extensions to the name of the R graphics device function that handles them,
# for extensions that do not literally match a function of the same name (e.g. "eps" is
# handled by grDevices::postscript(), not by a function called "eps"). Extensions not
# listed here are looked up as-is (get(extension)), as before.
exportgraphics.devicemap <-
  c (eps = "postscript",
     ps  = "postscript",
     jpg = "jpeg")

#' Open a graphics device
#'
#' Starts the graphics device driver
#' @name exportgraphics
#' @param file A character string giving the name of the file.
#' @param type The type of graphics device. Deduced from the file extension by default:
#' \code{"eps"} and \code{"ps"} go to \code{\link[grDevices]{postscript}}, \code{"jpg"} to
#' \code{\link[grDevices]{jpeg}}, and any other extension is taken as the name of an R
#' function (\code{"pdf"}, \code{"png"}, ...).
#' @param export If given, explicitly overrides the global export toggle set by \code{\link{toggleexport}} for
#' this call only. By default, the global toggle is used, so existing code is unaffected.
#' @param ... Other parameters.
#' @seealso \code{\link{closegraphics}}, \code{\link{toggleexport}}, \code{\link[grDevices]{Devices}}
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' # Extensions that don't match an R function name directly are now handled:
#' exportgraphics ("export.eps")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' }
exportgraphics <-
  function (file, type = tail (strsplit (file, split = "\\.") [[1]], 1), export = fdm2id.globals$export, ...)
  {
    if (is.character (type))
    {
      ext = tolower (type)
      fname = unname (exportgraphics.devicemap [ext])
      if (is.na (fname))
        fname = ext
      type = get (fname)
    }
    if (export)
      type (file, ...)
  }

#' @rdname toggleexport
#' @export
exportgraphics.off <-
  function ()
  {
    toggleexport (FALSE)
  }

#' @rdname toggleexport
#' @export
exportgraphics.on <-
  function ()
  {
    toggleexport (TRUE)
  }

#' Rotation
#'
#' Rotation on two variables of a numeric dataset
#' @name rotation
#' @param d The dataset.
#' @param angle The angle of the rotation.
#' @param axis The axis.
#' @param range The range of the angle (360, 2*pi, 100, ...)
#' @return A rotated data matrix.
#' @export
#' @examples
#' d = data.parabol ()
#' d [, -3] = rotation (d [, -3], 45, range = 360)
#' plotdata (d [, -3], d [, 3])
rotation = function (d, angle, axis = 1:2, range = 2 * pi)
{
  theta = 2 * pi * angle / range
  rot = diag (ncol (d))
  rot [axis, axis] = matrix (c (cos (theta), sin (theta), -sin (theta), cos (theta)), ncol = 2)
  res = as.matrix (d) %*% rot
  return (res)
}

#' Running time
#'
#' Return the running time of a function
#' @name runningtime
#' @param FUN The function to be evaluated.
#' @param ... The parameters to be passes to function \code{FUN}.
#' @return The running time of function \code{FUN}.
#' @export
#' @seealso \code{\link[base]{difftime}}
#' @examples
#' sqrt (x = 1:100)
#' runningtime (sqrt, x = 1:100)
runningtime <-
  function (FUN, ...)
  {
    start = Sys.time ()
    FUN (...)
    end = Sys.time ()
    return (end - start)
  }

#' Splits a dataset into training set and test set
#'
#' This function splits a dataset into training set and test set. Return an object of class \code{\link{dataset-class}}.
#' @name splitdata
#' @param dataset The dataset to be split (\code{data.frame} or \code{matrix}).
#' @param target The column index (numeric) or column name (character) of the target variable
#' (class label or response variable).
#' @param size The size of the training set: either a number of observations, or a proportion
#' between 0 and 1.
#' @param seed A specified seed for random number generation.
#' @param stratify Whether the split preserves the proportions of the classes. It matters as
#' soon as they are imbalanced: a plain random split can leave a rare class out of one side
#' altogether. Ignored when the target is numeric.
#' @return An object of class \code{\link{dataset-class}}.
#' @export
#' @seealso \code{\link{dataset-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' str (d)
splitdata <-
  function (dataset, target, size = round (0.7 * nrow (dataset)), seed = NULL, stratify = TRUE)
  {
    # A column name is accepted here too, as in augmentation(): the two functions sit side by
    # side in any practical session.
    target = column.index (dataset, target, "splitdata")
    setseed (seed)
    if (size < 1)
      size = round (size * nrow (dataset))
    y = dataset [, target]
    if (stratify && stratifiable (y))
      s = stratified.sample (y, size)
    else
      s = sample (nrow (dataset), size)
    train.x = dataset [s, -target]
    train.y = dataset [s, target]
    test.x = dataset [-s, -target]
    test.y = dataset [-s, target]
    res = list (train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y)
    class (res) = "dataset"
    return (res)
  }

#' Toggle graphic exports
#'
#' Toggle graphic exports on and off
#' @name toggleexport
#' @aliases exportgraphics.off exportgraphics.on toggleexport.off toggleexport.on
#' @param export If \code{TRUE}, exports are activated, if \code{FALSE}, exports are deactivated. If \code{null}, switches on and off.
#' @seealso \code{\link{closegraphics}}, \code{\link{exportgraphics}}
#' @rdname toggleexport
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' toggleexport (FALSE)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' toggleexport (TRUE)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' }
toggleexport <-
  function (export = NULL)
  {
    if (is.null (export))
      export = !fdm2id.globals$export
    fdm2id.globals$export = export
  }

#' @rdname toggleexport
#' @export
toggleexport.off <-
  function ()
  {
    toggleexport (FALSE)
  }

#' @rdname toggleexport
#' @export
toggleexport.on <-
  function ()
  {
    toggleexport (TRUE)
  }

#' Check and clean class labels
#'
#' Internal helper shared by the classification methods that cannot cope with empty classes.
#' It coerces \code{labels} to a \code{factor} (which makes it work with character vectors,
#' for which \code{nlevels} returns 0), drops the levels that are not observed (warning about
#' them) and checks that at least two classes remain.
#' @name check.classes
#' @param labels Class labels (\code{vector} or \code{factor}).
#' @param method The name of the calling method, used in the messages (a character string).
#' @param min.classes The minimal number of (non-empty) classes required (default: 2).
#' @return \code{labels}, as a \code{factor} with no empty level.
#' @keywords internal
check.classes <-
  function (labels, method = "", min.classes = 2)
  {
    prefix = if (nchar (method) > 0) paste0 (method, ": ") else ""
    # Note: factor (f) on a factor already drops the unused levels, so the empty ones must be
    # looked for *before* any coercion -- otherwise the warning below could never fire.
    ll = if (is.factor (labels)) labels else factor (labels)
    empty = levels (ll) [table (ll) == 0]
    if (length (empty) > 0)
    {
      warning (prefix, "the following class(es) have no observation and have been dropped: ",
               paste (empty, collapse = ", "), ".")
      ll = droplevels (ll)
    }
    if (nlevels (ll) < min.classes)
      stop (prefix, "at least ", min.classes, " non-empty classes are needed, but the labels ",
            "contain only ", nlevels (ll), " (", paste (levels (ll), collapse = ", "), ").")
    return (ll)
  }

#' @keywords internal
# Is stratification meaningful for this target? Only for a categorical one -- there is nothing
# to keep the proportions of in a numeric target.
stratifiable <-
  function (y) is.factor (y) || is.character (y) || is.logical (y)

#' @keywords internal
# Assigns each observation to one of 'nfolds' folds so that every class is spread over the
# folds in (as close as possible to) the proportions it has in the whole sample.
#
# Each class is shuffled and dealt round-robin; the starting fold carries over from one class
# to the next, so the leftover observations of the different classes do not all pile up in
# fold 1.
stratified.folds <-
  function (y, nfolds)
  {
    y = factor (y)
    folds = integer (length (y))
    offset = 0
    for (l in levels (y))
    {
      idx = which (y == l)
      idx = idx [sample (length (idx))]
      folds [idx] = ((offset + seq_along (idx) - 1) %% nfolds) + 1
      offset = offset + length (idx)
    }
    return (folds)
  }

#' @keywords internal
# Draws 'size' observations, keeping the class proportions. The per-class quotas are the exact
# proportions rounded down, the remaining places going to the classes with the largest
# fractional parts (largest-remainder method), so the quotas add up to exactly 'size'.
#
# A class with at least two observations always keeps at least one on each side of the split:
# a training set missing a class entirely cannot produce a model that predicts it, and a test
# set missing one measures nothing about it.
stratified.sample <-
  function (y, size)
  {
    y = factor (y)
    n = length (y)
    counts = as.vector (table (y))
    exact = size * counts / n
    # Preferred bounds: a class with at least two observations keeps at least one on each side
    # of the split. They are dropped when 'size' makes them impossible to satisfy -- asking for
    # 149 of 150 observations cannot leave one of each of three classes in the test set -- since
    # returning the requested number of observations matters more than the preference.
    lower = ifelse (counts >= 2, 1, 0)
    upper = ifelse (counts >= 2, counts - 1, counts)
    if ((size < sum (lower)) || (size > sum (upper)))
    {
      lower = rep (0, length (counts))
      upper = counts
    }
    # Largest-remainder allocation within those bounds.
    take = pmin (pmax (floor (exact), lower), upper)
    while (sum (take) < size)
    {
      room = which (take < upper)
      take [room [which.max ((exact - take) [room])]] =
        take [room [which.max ((exact - take) [room])]] + 1
    }
    while (sum (take) > size)
    {
      room = which (take > lower)
      take [room [which.min ((exact - take) [room])]] =
        take [room [which.min ((exact - take) [room])]] - 1
    }
    res = unlist (lapply (seq_along (levels (y)), function (i)
    {
      idx = which (y == levels (y) [i])
      return (idx [sample (length (idx), take [i])])
    }))
    # Shuffled, so that the training set is not ordered by class.
    return (res [sample (length (res))])
  }

#' @keywords internal
# The distance from every observation to its k-th nearest neighbour, itself excluded.
#
# Written as sort() of every column and then row k + 1, this sorted n values per column to
# read one of them: O(n^2 log n) where O(n^2) suffices. sort.int (partial = ) stops as soon as
# the wanted order statistic is in place.
#
# The distance matrix itself stays with flexclust::dist2(). Replacing it by the BLAS identity
# ||a - b||^2 = ||a||^2 + ||b||^2 - 2 a.b looked like an obvious win and measured as a loss
# (0.37 s against 0.20 s on 2000 x 6): dist2() is already C, and the identity needs three
# extra passes over an n x n matrix in R to assemble and clamp the result.
kdistances <-
  function (d, k)
  {
    dis = flexclust::dist2 (d, d)
    want = min (k + 1, nrow (dis))
    return (apply (dis, 2, function (v) sort.int (v, partial = want) [want]))
  }

