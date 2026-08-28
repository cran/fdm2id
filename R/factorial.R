#' Factorial analysis results
#'
#' This class contains the result of a factorial analysis, as obtained by \code{\link{CA}},
#' \code{\link{MCA}} or \code{\link{PCA}}.
#'
#' Objects of this class are the objects returned by the corresponding \pkg{FactoMineR}
#' functions (\code{\link[FactoMineR]{CA}}, \code{\link[FactoMineR]{MCA}} and
#' \code{\link[FactoMineR]{PCA}}), with the extra class \code{factorial} prepended so that
#' \code{\link{plot.factorial}} can provide a uniform plotting interface. The second class
#' (\code{"ca"}, \code{"mca"} or \code{"pca"}) records which analysis was performed.
#' @name factorial-class
#' @seealso \code{\link{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link{plot.factorial}}
NULL

#' @keywords internal
# How many axes a correspondence analysis of this table can produce.
#
# A contingency table of I rows by J columns carries min (I, J) - 1 non-trivial dimensions: the
# first one is the margins themselves, which carry no association and are removed. Supplementary
# rows and columns take no part in the axes, so they do not count.
ca.ncp <-
  function (d, row.sup = NULL, col.sup = NULL, quanti.sup = NULL, quali.sup = NULL)
  {
    rows = nrow (d) - length (row.sup)
    cols = ncol (d) - length (unique (c (col.sup, quanti.sup, quali.sup)))
    return (min (rows, cols) - 1)
  }

#' @keywords internal
# How many axes a multiple correspondence analysis of this table can produce.
#
# An MCA is a correspondence analysis of the indicator matrix, and the columns of a variable
# with J_q modalities sum to one: they carry J_q - 1 independent directions. Q active variables
# therefore give sum (J_q - 1), that is J - Q where J is the total number of modalities -- never
# more than n - 1.
#
# The modalities have to be counted as *observed* among the active individuals, not as declared:
# a level nobody takes carries nothing. On the first 200 rows of 'titanic' the declared levels
# promise six axes and the analysis produces two.
mca.ncp <-
  function (d, ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL)
  {
    active = setdiff (seq_len (ncol (d)), unique (c (quanti.sup, quali.sup)))
    rows = setdiff (seq_len (nrow (d)), ind.sup)
    if ((length (active) == 0) || (length (rows) == 0))
      return (0)
    modalities = sapply (active, function (j) nlevels (droplevels (factor (d [rows, j]))))
    return (min (sum (modalities - 1), length (rows) - 1))
  }

#' @keywords internal
# How many axes a principal component analysis of this table can produce: one per active
# variable, and never more than n - 1.
pca.ncp <-
  function (d, ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL)
  {
    variables = ncol (d) - length (unique (c (quanti.sup, quali.sup)))
    rows = nrow (d) - length (ind.sup)
    return (min (variables, rows - 1))
  }

#' @keywords internal
# Says so when there is no axis to compute, instead of letting FactoMineR fail inside its
# singular value decomposition on "max(nu, nv) must be positive".
factorial.checkncp <-
  function (ncp, fname, reason)
  {
    if (ncp < 1)
      stop (fname, ": there is no axis to compute -- ", reason, ".")
    return (ncp)
  }

#' Correspondence Analysis (CA)
#'
#' Performs Correspondence Analysis (CA) including supplementary row and/or column points.
#' @name CA
#' @param d A data frame or a table with n rows and p columns, i.e. a contingency table.
#' @param ncp The number of dimensions kept in the results. All of them, by default: a table of
#' I active rows by J active columns carries \code{min (I, J) - 1} of them, the first dimension
#' of a contingency table being its margins, which hold no association.
#' @param row.sup A vector indicating the indexes of the supplementary rows.
#' @param col.sup A vector indicating the indexes of the supplementary columns.
#' @param quanti.sup A vector indicating the indexes of the supplementary continuous variables.
#' @param quali.sup A vector indicating the indexes of the categorical supplementary variables.
#' @param row.w An optional row weights (by default, a vector of 1 for uniform row weights); the weights are given only for the active individuals.
#' @return The CA on the dataset.
#' @export
#' @seealso \code{\link[FactoMineR]{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link{plot.factorial}}, \code{\link{factorial-class}}
#' @examples
#' data (children, package = "FactoMineR")
#' CA (children, row.sup = 15:18, col.sup = 6:8)
CA <-
  function (d, ncp = ca.ncp (d, row.sup, col.sup, quanti.sup, quali.sup),
            row.sup = NULL, col.sup = NULL,
            quanti.sup = NULL, quali.sup = NULL, row.w = NULL)
  {
    ncp = factorial.checkncp (ncp, "CA",
                              "a table with a single active row or column holds no association")
    ca = FactoMineR::CA (d, ncp = ncp, row.sup = row.sup, col.sup = col.sup, quanti.sup = quanti.sup,
                         quali.sup = quali.sup, row.w = row.w, graph = FALSE)
    # Everything predict.factorial() needs to project new rows into this very analysis.
    ca$fdm2id = list (fun = FactoMineR::CA, data = d, supplementary = "row.sup",
                      args = list (ncp = ncp, row.sup = row.sup, col.sup = col.sup,
                                   quanti.sup = quanti.sup, quali.sup = quali.sup,
                                   row.w = row.w, graph = FALSE))
    class (ca) = c ("factorial", "ca", class (ca) [-1])
    return (ca)
  }

#' Kaiser rule
#'
#' Apply the Kaiser rule to determine the appropriate number of PCA axes.
#' @name kaiser
#' @param pca The PCA result (object of class \code{factorial-class}).
#' @export
#' @seealso \code{\link{PCA}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' pca = PCA (iris, quali.sup = 5)
#' kaiser (pca)
kaiser <-
  function (pca)
  {
    # The eigenvalues come out in decreasing order, so counting those above the average and
    # taking the rank of the last of them are the same thing -- except when there is none,
    # where max (which (...)) is max (integer (0)), i.e. -Inf with a warning. That happens on
    # uncorrelated standardised variables, whose correlation matrix is the identity: every
    # axis then carries exactly the average, and the rule keeps none.
    eig = pca$eig [, 1]
    return (sum (eig > mean (eig)))
  }

#' Multiple Correspondence Analysis (MCA)
#'
#' Performs Multiple Correspondence Analysis (MCA) with supplementary individuals, supplementary quantitative variables and supplementary categorical variables.
#' Performs also Specific Multiple Correspondence Analysis with supplementary categories and supplementary categorical variables.
#' Missing values are treated as an additional level, categories which are rare can be ventilated.
#' @name MCA
#' @param d A data frame or a table with n rows and p columns, i.e. a contingency table.
#' @param ncp The number of dimensions kept in the results. All of them, by default: one per
#' active variable, and never more than \eqn{n - 1}. A variable
#' with \eqn{J_q} modalities contributes \eqn{J_q - 1} of them -- its indicator columns sum to
#' one, so one of them is redundant -- which makes \eqn{J - Q} for \eqn{Q} active variables
#' holding \eqn{J} modalities between them, and never more than \eqn{n - 1}. The modalities
#' are counted as \emph{observed} among the active individuals: a level nobody takes carries
#' nothing.
#' @param ind.sup A vector indicating the indexes of the supplementary individuals.
#' @param quanti.sup A vector indicating the indexes of the quantitative supplementary variables.
#' @param quali.sup A vector indicating the indexes of the categorical supplementary variables.
#' @param row.w An optional row weights (by default, a vector of 1 for uniform row weights); the weights are given only for the active individuals.
#' @return The MCA on the dataset.
#' @export
#' @seealso \code{\link[FactoMineR]{MCA}}, \code{\link{CA}}, \code{\link{PCA}}, \code{\link{plot.factorial}}, \code{\link{factorial-class}}
#' @examples
#' data (tea, package = "FactoMineR")
#' MCA (tea, quanti.sup = 19, quali.sup = 20:36)
MCA <-
  function (d, ncp = mca.ncp (d, ind.sup, quanti.sup, quali.sup), ind.sup = NULL,
            quanti.sup = NULL, quali.sup = NULL, row.w = NULL)
  {
    ncp = factorial.checkncp (ncp, "MCA",
                              "every active variable takes a single value on these individuals")
    mca = FactoMineR::MCA (X = d, ncp = ncp, ind.sup = ind.sup, quanti.sup = quanti.sup,
                           quali.sup = quali.sup, row.w = row.w, graph = FALSE)
    mca$fdm2id = list (fun = FactoMineR::MCA, data = d, supplementary = "ind.sup",
                       args = list (ncp = ncp, ind.sup = ind.sup, quanti.sup = quanti.sup,
                                    quali.sup = quali.sup, row.w = row.w, graph = FALSE))
    class (mca) = c ("factorial", "mca", class (mca) [-1])
    return (mca)
  }

#' Principal Component Analysis (PCA)
#'
#' Performs Principal Component Analysis (PCA) with supplementary individuals, supplementary quantitative variables and supplementary categorical variables.
#' Missing values are replaced by the column mean.
#' @name PCA
#' @param d A data frame with n rows (individuals) and p columns (numeric variables).
#' @param scale.unit A boolean, if TRUE (value set by default) then data are scaled to unit variance.
#' @param ncp The number of dimensions kept in the results. All of them, by default: one per
#' active variable, and never more than \eqn{n - 1}.
#' @param ind.sup A vector indicating the indexes of the supplementary individuals.
#' @param quanti.sup A vector indicating the indexes of the quantitative supplementary variables.
#' @param quali.sup A vector indicating the indexes of the categorical supplementary variables.
#' @param row.w An optional row weights (by default, a vector of 1 for uniform row weights); the weights are given only for the active individuals.
#' @param col.w An optional column weights (by default, uniform column weights); the weights are given only for the active variables.
#' @return The PCA on the dataset.
#' @export
#' @seealso \code{\link[FactoMineR]{PCA}}, \code{\link{CA}}, \code{\link{MCA}}, \code{\link{plot.factorial}}, \code{\link{kaiser}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' PCA (iris, quali.sup = 5)
PCA <-
  function (d, scale.unit = TRUE, ncp = pca.ncp (d, ind.sup, quanti.sup, quali.sup),
            ind.sup = NULL,
            quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
            col.w = NULL)
  {
    ncp = factorial.checkncp (ncp, "PCA", "a single observation carries no variance")
    pca = FactoMineR::PCA (d, scale.unit = scale.unit [1], ncp = ncp, ind.sup = ind.sup, quanti.sup = quanti.sup,
                           quali.sup = quali.sup, row.w = row.w, col.w = col.w, graph = FALSE)
    pca$fdm2id = list (fun = FactoMineR::PCA, data = d, supplementary = "ind.sup",
                       args = list (scale.unit = scale.unit [1], ncp = ncp, ind.sup = ind.sup,
                                    quanti.sup = quanti.sup, quali.sup = quali.sup,
                                    row.w = row.w, col.w = col.w, graph = FALSE))
    class (pca) = c ("factorial", "pca", class (pca) [-1])
    return (pca)
  }

#' @keywords internal
plotfactorial.ind <-
  function (x, axes = c (1, 2), col = NULL, pch = NULL, labels = FALSE, legendpos = "topleft", ...)
  {
    coord = x$ind$coord [, axes]
    xlab = paste (colnames (coord) [1], " (", round (x$eig [axes [1], 2], 2), " %)", sep = "")
    ylab = paste (colnames (coord) [2], " (", round (x$eig [axes [2], 2], 2), " %)", sep = "")
    k = NULL
    if (is.null (col) && is.null (pch) && !is.null (x$call$quali.sup))
    {
      k = factor (x$call$quali.sup$quali.sup [, 1])
      col = as.numeric (k) + 1
      pch = as.numeric (k) + 1
    }
    if (is.null (col))
      col = 1
    if (is.null (pch))
      pch = 1
    if (labels)
    {
      graphics::plot (coord, col = 0, asp = 1, xlab = xlab, ylab = ylab, ...)
      graphics::text (coord, rownames (coord), col = col)
    }
    else
      graphics::plot (coord, col = col, pch = pch, asp = 1, xlab = xlab, ylab = ylab, ...)
    if (!is.null (x$ind.sup$coord))
    {
      coord.sup = x$ind.sup$coord [, axes]
      if (labels)
        graphics::text (coord.sup, rownames (coord.sup), col = "blue", font = 3)
      else
        graphics::points (coord.sup, col = "blue", pch = 3)
    }
    if (!is.null (k))
      graphics::legend (x = legendpos, legend = levels (k), pch = sort (unique (pch)), col = sort (unique (col)), bty = "n")
  }

#' Plot function for factorial-class
#'
#' Plot PCA, CA or MCA.
#' @name plot.factorial
#' @param x The PCA, CA or MCA result (object of class \code{factorial-class}).
#' @param type The graph to plot.
#' @param axes The factorial axes to be printed (numeric \code{vector}).
#' @param col Color(s) of the individuals (\code{type = "ind"}). If \code{NULL} (the default)
#' and a qualitative supplementary variable was given, the individuals are colored by it, as
#' \code{\link{plotdata}} does; otherwise any value the base \code{col} parameter accepts, so
#' that an external clustering can be used.
#' @param pch Point style(s) of the individuals on the scatter plot (\code{type = "ind"}, PCA only -- \code{\link[FactoMineR]{plot.CA}} and \code{\link[FactoMineR]{plot.MCA}}, which draw the CA and MCA plots, have no equivalent argument). Same default/auto-detection logic as \code{col}. Accepts the same values as the base \code{pch} graphical parameter.
#' @param labels Whether the row names are shown instead of points (\code{type = "ind"}).
#' @param legendpos Position of the legend (\code{type = "ind"}, PCA only, when individuals are colored by a qualitative variable).
#' @param ... Other parameters.
#' @method plot factorial
#' @export
#' @seealso \code{\link{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link{plotdata}}, \code{\link[FactoMineR]{plot.CA}}, \code{\link[FactoMineR]{plot.MCA}}, \code{\link[FactoMineR]{plot.PCA}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' pca = PCA (iris, quali.sup = 5)
#' plot (pca) # Automatically colored/legended by the qualitative supplementary variable
#' plot (pca, type = "cor")
#' plot (pca, type = "eig")
#' # Overriding colors and point styles manually (e.g. by an external clustering)
#' km = KMEANS (iris [, -5], k = 3)
#' plot (pca, col = km$cluster + 1, pch = km$cluster + 1)
plot.factorial <-
  function (x, type = c ("ind", "cor", "eig"), axes = c (1, 2), col = NULL, pch = NULL, labels = FALSE, legendpos = "topleft", ...)
  {
    if ("pca" %in% class (x))
      class (x) = c ("PCA", class (x) [-1])
    else if ("ca" %in% class (x))
      class (x) = c ("CA", class (x) [-1])
    else if ("mca" %in% class (x))
      class (x) = c ("MCA", class (x) [-1])
    if (type [1] == "ind")
    {
      if ("PCA" %in% class (x))
        plotfactorial.ind (x, axes = axes, col = col, pch = pch, labels = labels, legendpos = legendpos, ...)
      else if ("CA" %in% class (x))
      {
        # FactoMineR draws these plots itself, so 'col' and 'labels' are forwarded to the
        # arguments it understands rather than reimplemented ('pch' has no equivalent there).
        args = list (x = x, axes = axes, label = if (labels) "all" else "none")
        if (!is.null (col))
          args$col.row = col
        do.call (FactoMineR::plot.CA, c (args, list (...)))
      }
      else if ("MCA" %in% class (x))
      {
        args = list (x = x, choix = "ind", axes = axes, label = if (labels) "all" else "none")
        if (!is.null (col))
          args$col.ind = col
        do.call (FactoMineR::plot.MCA, c (args, list (...)))
      }
    }
    else if (type [1] == "cor")
    {
      if ("PCA" %in% class (x))
      {
        if (x$call$scale.unit)
          FactoMineR::plot.PCA (x, choix = "var", axes = axes, ...)
        else
          FactoMineR::plot.PCA (x, choix = "varcor", axes = axes, ...)
      }
      else if ("MCA" %in% class (x))
        FactoMineR::plot.MCA (x, choix = "var", axes = axes, ...)
      else
      {
        message ("Unavailable plot")
      }
    }
    else if (type [1] == "eig")
    {
      graphics::plot (x$eig [, 3], ylim = c (0, 100), col = 0, t = "b", xaxt = 'n', yaxt ='n',
                      xlab = "Axes", ylab = "Contribution", lwd = 2)
      graphics::grid ()
      graphics::lines (x = 1:nrow (x$eig), y = x$eig [, 3], type = "b", col = "red", lwd = 2)
      graphics::lines (x = 1:nrow (x$eig), y = x$eig [, 2], type = "b", col = "blue", lwd = 2)
      graphics::legend ("right", lty = 1, lwd = 2, col = c ("blue", "red"), bty = "n",
                        legend = c ("Percentage of variance", "Cumulative percentage of variance"))
      graphics::axis (side = 1, at = 1:nrow (x$eig), lwd = 0, lwd.ticks = 1)
      graphics::axis (side = 2, at = seq (0, to = 100, by = 20), lwd = 0, lwd.ticks = 1)
    }
    else
    {
      message ("Unavailable plot")
    }
  }

#' Plot function for factorial-class
#'
#' Print PCA, CA or MCA.
#' @name print.factorial
#' @param x The PCA, CA or MCA result (object of class \code{factorial-class}).
#' @param ... Other parameters.
#' @method print factorial
#' @export
#' @seealso \code{\link{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link[FactoMineR]{print.CA}}, \code{\link[FactoMineR]{print.MCA}}, \code{\link[FactoMineR]{print.PCA}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' pca = PCA (iris, quali.sup = 5)
#' print (pca)
print.factorial <-
  function (x, ...)
  {
    if ("pca" %in% class (x))
    {
      class (x) = c ("PCA", class (x) [-1])
      FactoMineR::print.PCA (x, ...)
    }
    else if ("ca" %in% class (x))
    {
      class (x) = c ("CA", class (x) [-1])
      FactoMineR::print.CA (x, ...)
    }
    else if ("mca" %in% class (x))
    {
      class (x) = c ("MCA", class (x) [-1])
      FactoMineR::print.MCA (x, ...)
    }
  }

#' Projection of new observations into a factorial space
#'
#' Projects new observations into the factorial space computed by \code{\link{CA}},
#' \code{\link{MCA}} or \code{\link{PCA}} -- the same operation that is applied to the
#' observations the analysis was fitted on: centering (and, for \code{\link{PCA}} with
#' \code{scale.unit = TRUE}, scaling) with the parameters \emph{of the training data}, then
#' projection on the axes already computed. The axes are not recomputed and the new
#' observations have no influence on them: they are \emph{supplementary} individuals.
#'
#' The projection is obtained by handing the new rows back to \pkg{FactoMineR} as supplementary
#' individuals of the original analysis, so the coordinates are exactly the ones
#' \code{PCA (rbind (train, test), ind.sup = ...)} would give. The active analysis is refitted
#' in the process, which is unnoticeable on the sizes this package is meant for.
#'
#' Supplementary variables (\code{quanti.sup}, \code{quali.sup}) play no part in the axes, so
#' \code{test} does not have to carry them: any column of the training data that is missing from
#' \code{test} is filled in (with the training mean, or the first level) purely so that the two
#' can be stacked.
#' @name predict.factorial
#' @param object The factorial analysis (object of class \code{\link{factorial-class}}).
#' @param test The new observations, a \code{data.frame} or \code{matrix} with the same
#' (active) variables as the data the analysis was fitted on.
#' @param ... Other parameters.
#' @return The coordinates of the new observations on the factorial axes (a \code{matrix},
#' one row per observation and one column per axis).
#' @export
#' @method predict factorial
#' @seealso \code{\link{PCA}}, \code{\link{CA}}, \code{\link{MCA}},
#' \code{\link{factorial-class}}, \code{\link{predict.cda}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' pca = PCA (d$train.x)
#' # The coordinates of unseen observations on the axes of the training analysis
#' head (predict (pca, d$test.x))
#' # An observation of the training set projects onto the coordinates the analysis gave it
#' pca$ind$coord [1, ]
#' predict (pca, d$train.x [1, ])
predict.factorial <-
  function (object, test, ...)
  {
    kept = object$fdm2id
    if (is.null (kept))
      stop ("predict.factorial: this object was not built by fdm2id's CA(), MCA() or PCA() ",
            "-- it does not carry the information needed to project new observations.")
    train = kept$data
    test = as.data.frame (test)
    if (is.null (colnames (test)) && (ncol (test) == ncol (train)))
      colnames (test) = colnames (train)
    unknown = setdiff (colnames (test), colnames (train))
    if (length (unknown) > 0)
      stop ("predict.factorial: 'test' has ", length (unknown), " variable(s) the analysis ",
            "does not know: ", paste (unknown, collapse = ", "), ".")
    # Columns absent from 'test' can only be supplementary ones, which take no part in the
    # axes; they are filled with a harmless value so that the two tables can be stacked.
    for (v in setdiff (colnames (train), colnames (test)))
      test [[v]] = if (is.numeric (train [[v]])) mean (train [[v]], na.rm = TRUE)
                   else factor (levels (factor (train [[v]])) [1],
                                levels = levels (factor (train [[v]])))
    test = test [, colnames (train), drop = FALSE]
    for (v in colnames (train))
      if (is.factor (train [[v]]))
        test [[v]] = factor (as.character (test [[v]]), levels = levels (train [[v]]))
    both = rbind (train, test)
    rownames (both) = make.unique (as.character (c (rownames (train), rownames (test))))
    new = nrow (train) + seq_len (nrow (test))
    args = kept$args
    # An analysis fitted with supplementary rows of its own keeps them supplementary.
    args [[kept$supplementary]] = sort (unique (c (args [[kept$supplementary]], new)))
    res = do.call (kept$fun, c (list (both), args))
    coord = if (kept$supplementary == "row.sup") res$row.sup$coord else res$ind.sup$coord
    keep = match (as.character (rownames (both) [new]), rownames (coord))
    return (coord [keep, , drop = FALSE])
  }
