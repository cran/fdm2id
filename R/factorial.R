#' Factorial analysis results
#'
#' This class contains the classification model obtained by the CDA method.
#' @name factorial-class
#' @exportClass factorial
#' @seealso \code{\link{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link{plot.factorial}}
setClass ("factorial", representation ())

#' Correspondence Analysis (CA)
#'
#' Performs Correspondence Analysis (CA) including supplementary row and/or column points.
#' @name CA
#' @param d A ddata frame or a table with n rows and p columns, i.e. a contingency table.
#' @param ncp The number of dimensions kept in the results (by default 5).
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
  function (d, ncp = 5, row.sup = NULL, col.sup = NULL,
            quanti.sup = NULL, quali.sup = NULL, row.w = NULL)
  {
    ca = FactoMineR::CA (d, ncp = ncp, row.sup = row.sup, col.sup = col.sup, quanti.sup = quanti.sup,
                           quali.sup = quali.sup, row.w = row.w, graph = FALSE)
    class (ca) = c ("factorial", "ca", class (ca) [-1])
    return (ca)
  }

#' Keiser rule
#'
#' Apply the keiser rule to determine the appropriate number of PCA axes.
#' @name keiser
#' @param pca The PCA result (object of class \code{factorial-class}).
#' @export
#' @seealso \code{\link{PCA}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' pca = PCA (iris, quali.sup = 5)
#' keiser (pca)
keiser <-
  function (pca)
  {
    return (max (which (pca$eig [, 1] > mean (pca$eig [, 1]))))
  }

#' Multiple Correspondence Analysis (MCA)
#'
#' Performs Multiple Correspondence Analysis (MCA) with supplementary individuals, supplementary quantitative variables and supplementary categorical variables.
#' Performs also Specific Multiple Correspondence Analysis with supplementary categories and supplementary categorical variables.
#' Missing values are treated as an additional level, categories which are rare can be ventilated.
#' @name MCA
#' @param d A ddata frame or a table with n rows and p columns, i.e. a contingency table.
#' @param ncp The number of dimensions kept in the results (by default 5).
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
  function (d, ncp = 5, ind.sup = NULL,
            quanti.sup = NULL, quali.sup = NULL, row.w = NULL)
  {
    mca = FactoMineR::MCA (X = d, ncp = ncp, ind.sup = ind.sup, quanti.sup = quanti.sup,
                         quali.sup = quali.sup, row.w = row.w, graph = FALSE)
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
#' @param ncp The number of dimensions kept in the results (by default 5).
#' @param ind.sup A vector indicating the indexes of the supplementary individuals.
#' @param quanti.sup A vector indicating the indexes of the quantitative supplementary variables.
#' @param quali.sup A vector indicating the indexes of the categorical supplementary variables.
#' @param row.w An optional row weights (by default, a vector of 1 for uniform row weights); the weights are given only for the active individuals.
#' @param col.w An optional column weights (by default, uniform column weights); the weights are given only for the active variables.
#' @return The PCA on the dataset.
#' @export
#' @seealso \code{\link[FactoMineR]{PCA}}, \code{\link{CA}}, \code{\link{MCA}}, \code{\link{plot.factorial}}, \code{\link{keiser}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' PCA (iris, quali.sup = 5)
PCA <-
  function (d, scale.unit = TRUE, ncp = 5, ind.sup = NULL,
                 quanti.sup = NULL, quali.sup = NULL, row.w = NULL,
                 col.w = NULL)
  {
    pca = FactoMineR::PCA (d, scale.unit = scale.unit, ncp = ncp, ind.sup = ind.sup, quanti.sup = quanti.sup,
                     quali.sup = quali.sup, row.w = row.w, col.w = col.w, graph = FALSE)
    class (pca) = c ("factorial", "pca", class (pca) [-1])
    return (pca)
  }

#' Plot function for factorial-class
#'
#' Plot PCA, CA or MCA.
#' @name plot.factorial
#' @param x The PCA, CA or MCA result (object of class \code{factorial-class}).
#' @param type The graph to plot.
#' @param axes The factorial axes to be printed (numeric \code{vector}).
#' @param ... Other parameters.
#' @method plot factorial
#' @export
#' @seealso \code{\link{CA}}, \code{\link{MCA}}, \code{\link{PCA}}, \code{\link[FactoMineR]{plot.CA}}, \code{\link[FactoMineR]{plot.MCA}}, \code{\link[FactoMineR]{plot.PCA}}, \code{\link{factorial-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' pca = PCA (iris, quali.sup = 5)
#' plot (pca)
#' plot (pca, type = "cor")
#' plot (pca, type = "eig")
plot.factorial <-
  function (x, type = c ("ind", "cor", "eig"), axes = c (1, 2), ...)
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
        FactoMineR::plot.PCA (x, choix = "ind", axes = axes, ...)
      else if ("CA" %in% class (x))
        FactoMineR::plot.CA (x, axes = axes, ...)
      else if ("MCA" %in% class (x))
        FactoMineR::plot.MCA (x, choix = "ind", axes = axes, ...)
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
