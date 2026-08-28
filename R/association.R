#' APRIORI classification model
#'
#' This class contains the classification model obtained by the APRIORI association rules method.
#'
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{rules}}{The set of rules obtained by APRIORI.}
#'   \item{\code{transactions}}{The training set as a \code{transaction} object.}
#'   \item{\code{train}}{The training set (description). A \code{matrix} or \code{data.frame}.}
#'   \item{\code{labels}}{Class labels of the training set. Either a \code{factor} or an integer \code{vector}.}
#'   \item{\code{supp}}{The minimal support of an item set (numeric value).}
#'   \item{\code{conf}}{The minimal confidence of an item set (numeric value).}
#' }
#' @name apriori-class
#' @seealso \code{\link{APRIORI}}, \code{\link{predict.apriori}}, \code{\link{print.apriori}},
#' \code{\link{summary.apriori}}, \code{\link[arules]{apriori}}
NULL

#' Classification using APRIORI
#'
#' This function builds a classification model using the association rules method APRIORI.
#' @name APRIORI
#' @param train The training set (description), as a \code{data.frame}.
#' @param labels Class labels of the training set (\code{vector} or \code{factor}).
#' @param supp The minimal support of an item set (numeric value).
#' @param conf The minimal confidence of an item set (numeric value).
#' @param prune A logical indicating whether to prune redundant rules or not (default: \code{FALSE}).
#' @inheritParams tune.doc
#' @param methodparameters Present for interface consistency with \code{\link{performance}}
#' (which always passes it when fitting a model). Currently unused: \code{APRIORI} does not
#' support reusing pre-tuned parameters.
#' @param graph Present for interface consistency with \code{\link{performance}} (which always
#' passes it when fitting a model). Currently unused: \code{APRIORI} does not produce a plot.
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
  function (train, labels, supp = .05, conf = .8, prune = FALSE,
            tune = FALSE, methodparameters = NULL, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    res = NULL
    if (tune)
      res = emptyparams ()
    else
    {
      ls = cbind.data.frame (train, Class = labels)
      tr = methods::as (ls, "transactions")
      # '...' was declared and went nowhere, so arules' own parameters -- 'maxlen' above all,
      # whose default of 10 truncates the search and warns about it -- could not be reached.
      parameter = list (supp = supp, conf = conf, minlen = 1)
      extra = list (...)
      parameter [names (extra)] = extra
      apr = arules::apriori (tr, parameter = parameter, control = list (verbose = FALSE))
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
#' r = apriori (Adult, parameter = list (supp = .4, conf = .8))
#' inspect (filter.rules (r, right = "marital-status="))
#' # The equivalent call in arules itself
#' subset (r, subset = rhs %pin% "marital-status=")
filter.rules <-
  function (rules, pattern = NULL, left = pattern, right = pattern, removeMatches = FALSE)
  {
    if (is.null (left) && is.null (right))
      return (rules)
    matched = rep (FALSE, length (rules))
    if (!is.null (left))
      matched = matched | arules::`%pin%` (arules::lhs (rules), left)
    if (!is.null (right))
      matched = matched | arules::`%pin%` (arules::rhs (rules), right)
    if (removeMatches)
      return (rules [!matched])
    else
      return (rules [matched])
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
#' # The default support (0.1) yields ~6000 rules on Adult, and general.rules() compares every
#' # pair of them twice: that single call took more than 8 seconds. A higher support keeps the
#' # example instructive (169 rules, of which 8 are general) and instantaneous.
#' r = apriori (Adult, parameter = list (supp = .4, conf = .8))
#' inspect (general.rules (r))
general.rules <-
  function (r)
  {
    # arules::is.subset() answers with a *sparse* matrix; as.matrix() used to expand it, which
    # for the 6137 rules the default support yields on 'Adult' is a 37-million-cell logical
    # matrix built to hold a few thousand TRUE. Kept sparse, and with the diagonal subtracted
    # from the column sums rather than overwritten, the same 18 rules come out in half the time.
    subsets = arules::is.subset (r@lhs, r@lhs) & arules::is.subset (r@rhs, r@rhs)
    redundant = (Matrix::colSums (subsets) - as.numeric (Matrix::diag (subsets))) >= 1
    return (r [!redundant])
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
    supp = object$rules@quality$support
    rhs = factor (labels (r@rhs, setStart = "", setEnd = ""))
    n = nlevels (rhs) + 1
    l = gsub ("Class=", "", levels (rhs), fixed = TRUE)
    if (unmatched %in% l)
      stop ("predict.apriori: the 'unmatched' label (\"", unmatched, "\") collides with an ",
            "existing class label. Please choose another value for 'unmatched'.")
    rhs = as.numeric (rhs)
    lhs = methods::as (r@lhs, "matrix")
    lhs = lhs [, colnames (methods::as (t, "matrix"))]
    lhs = methods::as (lhs, "itemMatrix")
    pred = apply (arules::is.subset (lhs, t), 2, select.rule, rhs, q, supp, n)
    # select.rule() codes an unmatched observation as n, one past the last consequent, so 'l'
    # must always have n entries for l [code] to mean anything.
    l = c (l, unmatched)
    codes = sort (unique (pred))
    pred = factor (pred, levels = codes, labels = l [codes])
    return (pred)
  }

#' Plot function for apriori-class
#'
#' Plot the association rules obtained by APRIORI, using arulesViz.
#' @name plot.apriori
#' @param x The classification model (object of class \code{apriori-class}, created by \code{APRIORI}).
#' @param method The type of plot (see \code{\link[arulesViz]{plot}}, e.g. "scatterplot", "graph", "grouped", "paracoord").
#' @param measure,shading Parameters passed to \code{\link[arulesViz]{plot}}.
#' @param ... Other parameters passed to \code{\link[arulesViz]{plot}}.
#' @method plot apriori
#' @export
#' @seealso \code{\link{APRIORI}}, \code{\link{apriori-class}}, \code{\link[arulesViz]{plot}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' d = discretizeDF (iris,
#'     default = list (method = "interval", breaks = 3, labels = c ("small", "medium", "large")))
#' model = APRIORI (d [, -5], d [, 5], supp = .1, conf = .9, prune = TRUE)
#' plot (model)
#' plot (model, method = "graph")
#' }
plot.apriori <-
  function (x, method = "scatterplot", measure = c ("support", "confidence"), shading = "lift", ...)
  {
    # arulesViz does not export a 'plot' function -- it only registers an S3 method
    # (plot.rules) on the base plot() generic (see its NAMESPACE: S3method(plot,rules), no
    # export(plot)). arulesViz::plot(...) is therefore invalid ("object not exported"); the
    # plain generic below dispatches correctly since arulesViz is in Depends.
    plot (x$rules, method = method, measure = measure, shading = shading, ...)
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

#' @keywords internal
# The consequent of the rule an observation is classified by: the most confident of the rules
# whose antecedent it matches, ties broken by support and then by the order apriori() produced
# the rules in. 'n' -- one past the last consequent -- codes an observation no rule matches.
select.rule <-
  function (v, r, conf, supp, n)
  {
    if (!any (v))
      return (n)
    l = which (v)
    l = l [conf [l] == max (conf [l])]
    l = l [supp [l] == max (supp [l])]
    return (r [l [1]])
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
