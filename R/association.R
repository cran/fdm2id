#' APRIORI classification model
#'
#' This class contains the classification model obtained by the APRIORI association rules method.
#'
#' @name apriori-class
#' @slot rules The set of rules obtained by APRIORI.
#' @slot transactions The training set as a \code{transaction} object.
#' @slot train The training set (description). A \code{matrix} or \code{data.frame}.
#' @slot labels Class labels of the training set. Either a \code{factor} or an integer \code{vector}.
#' @slot supp The minimal support of an item set (numeric value).
#' @slot conf The minimal confidence of an item set (numeric value).
#' @exportClass apriori
#' @seealso \code{\link{APRIORI}}, \code{\link{predict.apriori}}, \code{\link{print.apriori}},
#' \code{\link{summary.apriori}}, \code{\link[arules]{apriori}}
setClass ("apriori",
          representation (rules = "ANY",
                          transactions = "ANY",
                          train = "data.frame",
                          labels = "vector",
                          supp = "numeric",
                          conf = "numeric"))

#' Classification using APRIORI
#'
#' This function builds a classification model using the association rules method APRIORI.
#' @name APRIORI
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param supp The minimal support of an item set (numeric value).
#' @param conf The minimal confidence of an item set (numeric value).
#' @param prune A logical indicating whether to prune redundant rules or not (default: \code{FALSE}).
#' @param tune If true, the function returns paramters instead of a classification model.
#' @param ... Other parameters.
#' @return The classification model, as an object of class \code{apriori}.
#' @export
#' @seealso \code{\link{predict.apriori}}, \code{\link{apriori-class}}, \code{\link[arules]{apriori}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = discretizeDF (iris,
#'     default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
#' APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
APRIORI <-
  function (train, labels, supp = .05, conf = .8, prune = FALSE, tune = FALSE, ...)
  {
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      ls = cbind.data.frame (train, Class = labels)
      tr = methods::as (ls, "transactions")
      apr = arules::apriori (tr, parameter = list (supp = supp, conf = conf, minlen = 1), control = list (verbose = FALSE))
      apr = filter.rules (apr, right = "Class=")
      if (prune)
        apr = general.rules (apr)
      res = list (rules = apr, transactions = tr, train = train, labels = labels, supp = supp, conf = conf)
      class (res) = "apriori"
    }
    return (res)
  }

#' Filtering a set of rules
#'
#' This function facilitate the selection of a subset from a set of rules.
#' @name filter.rules
#' @param rules A set of rules.
#' @param pattern A pattern to match (antecedent and consequent): a character string.
#' @param left A pattern to match (antecedent only): a character string.
#' @param right A pattern to match (consequent only): a character string.
#' @param removeMatches A logical indicating whether to remove matching rules (\code{TRUE}) or to keep those (\code{FALSE}).
#' @return The filtered set of rules.
#' @export
#' @seealso \code{\link[arules]{apriori}}, \code{\link[arules]{subset}}
#' @examples
#' require ("arules")
#' data ("Adult")
#' r = apriori (Adult)
#' filter.rules (r, right = "marital-status=")
#' subset (r, subset = rhs %pin% "marital-status=")
filter.rules <-
  function (rules, pattern = NULL, left = pattern, right = pattern, removeMatches = FALSE)
  {
    res = rules
    if (removeMatches)
    {
      tmp = NULL
      if (!is.null (left))
        res = arules::subset (rules, subset = arules::`%pin%` (lhs, left))
      if (!is.null (right))
      {
        tmp = arules::subset (rules, subset = arules::`%pin%` (rhs, right))
        if (!is.null (left))
          res = c (res, tmp)
        else
          res = tmp
      }
    }
    else
    {
      tmp = NULL
      if (!is.null (left))
        res = arules::subset (rules, subset = arules::`%pin%` (lhs, left))
      if (!is.null (right))
      {
        tmp = arules::subset (rules, subset = arules::`%pin%` (rhs, right))
        if (!is.null (left))
          res = c (res, tmp)
        else
          res = tmp
      }
    }
    return (res)
  }

#' Remove redundancy in a set of rules
#'
#' This function remove every redundant rules, keeping only the most general ones.
#' @name general.rules
#' @param r A set of rules.
#' @return A set of rules, without redundancy.
#' @export
#' @seealso \code{\link[arules]{apriori}}
#' @examples
#' require ("arules")
#' data ("Adult")
#' r = apriori (Adult)
#' inspect (general.rules (r))
general.rules <-
  function (r)
  {
    subset.matrix <- as.matrix (arules::is.subset(r@lhs, r@lhs) & arules::is.subset (r@rhs, r@rhs))
    diag (subset.matrix) = FALSE
    redundant <- colSums (subset.matrix, na.rm = TRUE) >= 1
    r <- r [!redundant]
    return (r)
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained by \code{apriori.classif}.
#' Observations that do not match any of the rules are labelled as "unmatched".
#' @name predict.apriori
#' @param object The classification model (of class \code{apriori}, created by \code{apriori.classif}).
#' @param test The test set (a \code{data.frame})
#' @param unmatched The class label given to the unmatched observations (a character string).
#' @param ... Other parameters.
#' @return A vector of predicted values (\code{factor}).
#' @export
#' @method predict apriori
#' @seealso \code{\link{APRIORI}}, \code{\link{apriori-class}}, \code{\link[arules]{apriori}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = discretizeDF (iris,
#'     default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
#' model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
#' predict (model, d [, -5])
predict.apriori <-
  function (object, test, unmatched  = "Unknown", ...)
  {
    t = methods::as (test, "transactions")
    r = object$rules
    q = object$rules@quality$confidence
    rhs = factor (labels (r@rhs, setStart = "", setEnd = ""))
    n = nlevels (rhs) + 1
    l = gsub ("Class=", "", levels (rhs), fixed = TRUE)
    rhs = as.numeric (rhs)
    lhs = methods::as (r@lhs, "matrix")
    lhs = lhs [, colnames (methods::as (t, "matrix"))]
    lhs = methods::as (lhs, "itemMatrix")
    pred = apply (arules::is.subset (lhs, t), 2, select.rule, rhs, q, n)
    if (length (unique (pred)) == n)
      l = c (l, unmatched)
    l = l [sort (unique (pred))]
    pred = factor (pred, labels = l)
    return (pred)
  }

#' Print a classification model obtained by APRIORI
#'
#' Print the set of rules in the classification model.
#' @name print.apriori
#' @param x The model to be printed.
#' @param ... Other parameters.
#' @export
#' @method print apriori
#' @seealso \code{\link{APRIORI}}, \code{\link{predict.apriori}}, \code{\link{summary.apriori}},
#' \code{\link{apriori-class}}, \code{\link[arules]{apriori}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = discretizeDF (iris,
#'     default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
#' model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
#' print (model)
print.apriori <-
  function (x, ...) arules::inspect (x$rules, ...)

select.rule <-
  function (v, r, q, n)
  {
    res = n
    if (any (v))
    {
      l = which (v)
      pred = r [l]
      conf = q [l]
      s = which (conf == max (conf))
      res = pred [s] [sample (length (s), 1)]
    }
    return (res)
  }

#' Print summary of a classification model obtained by APRIORI
#'
#' Print summary of the set of rules in the classification model obtained by APRIORI.
#' @name summary.apriori
#' @param object The model to be printed.
#' @param ... Other parameters.
#' @export
#' @method summary apriori
#' @seealso \code{\link{APRIORI}}, \code{\link{predict.apriori}}, \code{\link{print.apriori}},
#' \code{\link{apriori-class}}, \code{\link[arules]{apriori}}
#' @examples
#' require ("datasets")
#' data (iris)
#' d = discretizeDF (iris,
#'     default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
#' model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
#' summary (model)
summary.apriori <-
  function (object, ...) summary (object$rules, ...)
