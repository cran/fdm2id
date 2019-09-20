#' Non-negative Matrix Factorization
#'
#' Return the NMF decomposition.
#' @name NMF
#' @param x A numeric dataset (data.frame or matrix).
#' @param rank Specification of the factorization rank.
#' @param nstart How many random sets should be chosen?
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[NMF]{nmf}}
#' @examples
#' require (datasets)
#' data (iris)
#' NMF (iris [, -5], nstart = 1)
NMF <-
  function (x, rank = 2, nstart = 10, ...)
  {

    res = NMF::nmf (x, rank)
    eval = res@residuals
    for (i in 1:(nstart - 1))
    {
      tmp = NMF::nmf (x, rank)
      if (tmp@residuals < eval)
      {
        res = tmp
        eval = res@residuals
      }
    }
    colnames (res@fit@W) = paste ("Dim.", 1:rank)
    return (res)
  }

#' @keywords internal
panel.hist <- function(x, ...)
{
  usr = graphics::par ("usr")
  on.exit (graphics::par (usr))
  graphics::par (usr = c (usr [1:2], 0, 1.5))
  h = graphics::hist (x, plot = FALSE)
  breaks = h$breaks
  nB = length (breaks)
  y = h$counts
  y = y / max (y)
  graphics::rect (breaks [-nB], 0, breaks [-1], y, col="grey")
}

#' @keywords internal
panel.blank <- function(x, y) {}

#' @keywords internal
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr = graphics::par ("usr")
  on.exit (graphics::par (usr))
  graphics::par (usr = c (0, 1, 0, 1))
  r = stats::cor (x, y)
  txt = format (c (r, 0.123456789), digits=digits) [1]
  txt = paste(prefix, txt, sep = "")
  if (missing (cex.cor))
    cex = 0.6 / graphics::strwidth (txt)
  graphics::text (0.5, 0.5, txt, cex = cex)
}

#' Advanced plot function
#'
#' Plot a dataset.
#' @name plotdata
#' @param d A numeric dataset (data.frame or matrix).
#' @param k A categorial variable (vector or factor).
#' @param type The type of graphic to be plotted.
#' @param legendpos Position of the legend
#' @param alpha Color opacity (0-255).
#' @param asp Aspect ratio (default: 1).
#' @param labels Indicates whether or not labels (row names) should be showned on the (scatter) plot.
#' @param ... Other parameters.
#' @export
#' @examples
#' require (datasets)
#' data (iris)
#' # Without classification
#' plotdata (iris [, -5]) # Défault (pairs)
#' # With classification
#' plotdata (iris [, -5], iris [, 5]) # Défault (pairs)
#' plotdata (iris [, -5], iris [, 5], type = "scatter") # Scatter plot (PCA axis)
#' plotdata (iris [, -5], iris [, 5], type = "boxplot") # Boxplot
#' plotdata (iris [, -5], iris [, 5], type = "pca") # Scatter plot (PCA axis)
#' plotdata (iris [, -5], iris [, 5], type = "cda") # Scatter plot (CDA axis)
#' plotdata (iris [, -5], iris [, 5], type = "svd") # Scatter plot (SVD axis)
#' plotdata (iris [, -5], iris [, 5], type = "som") # Kohonen map
#' # With only one variable
#' plotdata (iris [, 1], iris [, 5]) # Défault (data vs. index)
#' plotdata (iris [, 1], iris [, 5], type = "scatter") # Scatter plot (data vs. index)
#' plotdata (iris [, 1], iris [, 5], type = "boxplot") # Boxplot
#' # With two variables
#' plotdata (iris [, 3:4], iris [, 5]) # Défault (scatter plot)
#' plotdata (iris [, 3:4], iris [, 5], type = "scatter") # Scatter plot
plotdata <-
  function (d, k = NULL, type = c ("pairs", "scatter", "boxplot", "pca", "cda", "svd", "nmf", "tsne", "som"), legendpos = "topleft", alpha = 200, asp = 1, labels = FALSE, ...)
  {
    col = 1
    pch = 1
    if (!is.null (k))
    {
      if (!is.factor (k))
        k = factor (k)
      pch = as.numeric (k)
      col = pch + 1
    }
    lcol = sort (unique (col))
    lpch = sort (unique (pch))
    if (is.vector (d))
    {
      if ((type [1] == "scatter") | (type [1] == "pairs"))
      {
        graphics::plot (cbind (Index = 1:(length (d)), Data = d), col = col, pch = pch, ...)
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "boxplot")
      {
        mini = min (d)
        maxi = max (d)
        names = colnames (d)
        v = 1
        if (!is.null (k))
        {
          graphics::boxplot (d~k, ylim = c (mini, maxi), ylab = "", col = lcol, xaxt='n', xlab = "")
          graphics::legend (x = legendpos, legend = levels (k), fill = lcol, bty = "n")
        }
        else
          graphics::boxplot (d, ylim = c (mini, maxi), ylab = "", col = "grey", xaxt='n', xlab = "")
      }
      else
        message ("Unavailable plot")
    }
    else
    {
      if ((type [1] == "pairs") & (ncol (d) > 2))
        graphics::pairs (d, upper.panel = panel.cor, diag.panel = panel.hist, pch = pch, col = col, asp = asp, ...)
      else if ((type [1] == "scatter") | (type [1] == "pairs"))
      {
        dd = NULL
        if (ncol (d) == 2)
          dd = d
        else
          dd = FactoMineR::PCA (d, scale.unit = FALSE, ncp = 2, graph = FALSE)$ind$coord [, 1:2]
        if (labels)
        {
          graphics::plot (dd, col = 0, asp = 1, ...)
          graphics::text (dd, row.names (d), col = col)
        }
        else
          graphics::plot (dd, col = col, pch = pch, asp = 1, ...)
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "boxplot")
      {
        mini = min (d)
        maxi = max (d)
        names = colnames (d)
        v = 1
        if (!is.null (k))
        {
          nbclusters = length (unique (col))
          nbclusters * ncol (d)
          at = (0:(ncol (d) - 1)) * nbclusters + 1 + (nbclusters - 1) / 2
          v = (1:(ncol (d) - 1)) * nbclusters + .5
          d = utils::stack (d)
          d$cluster = k
          graphics::boxplot (values~cluster+ind, d, ylim = c (mini, maxi), ylab = "", col = 2:(nbclusters + 1), xaxt='n', xlab = "")
        }
        else
        {
          graphics::boxplot (d, ylim = c (mini, maxi), ylab = "", col = "grey", xaxt='n', xlab = "")
          at = 1:ncol (d)
        }
        graphics::axis (side = 1, at = at, labels = names, lwd.ticks = FALSE, lwd = 0)
        if (!is.null (k))
        {
          graphics::abline (v = v, lty = 2, col = "grey")
          graphics::legend (x = legendpos, legend = levels (k), fill = lcol, bty = "n")
        }
      }
      else if (type [1] == "pca")
      {
        dd = FactoMineR::PCA (d, scale.unit = FALSE, ncp = 2, graph = FALSE)$ind$coord [, 1:2]
        if (labels)
        {
          graphics::plot (dd, col = 0, asp = 1, ...)
          graphics::text (dd, row.names (d), col = col)
        }
        else
          graphics::plot (dd, col = col, pch = pch, asp = 1, ...)
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "cda")
      {
        if (is.null (k))
          message ("Unavailable plot")
        else
        {
          dd = CDA (d, k)$proj [, 1:2]
          if (labels)
          {
            graphics::plot (dd, col = 0, asp = 1, ...)
            graphics::text (dd, row.names (d), col = col)
          }
          else
            graphics::plot (dd, col = col, pch = pch, asp = 1, ...)
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
        }
      }
      else if (type [1] == "svd")
      {
        dd = SVD (d)$proj$ind [, 1:2]
        if (labels)
        {
          graphics::plot (dd, col = 0, asp = 1, ...)
          graphics::text (dd, row.names (d), col = col)
        }
        else
          graphics::plot (dd, col = col, pch = pch, asp = 1, ...)
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "nmf")
      {
        dd = NMF (d)@fit@W [, 1:2]
        if (labels)
        {
          graphics::plot (dd, col = 0, asp = 1, ...)
          graphics::text (dd, row.names (d), col = col)
        }
        else
          graphics::plot (dd, col = col, pch = pch, asp = 1, ...)
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "tsne")
      {
        dd = TSNE (d, ...)$Y
        if (labels)
        {
          withCallingHandlers (suppressWarnings (graphics::plot (dd, col = 0, asp = 1, ...)), warning = function(w) {message(w)})
          graphics::text (dd, row.names (d), col = col)
        }
        else
          withCallingHandlers (suppressWarnings (graphics::plot (dd, col = col, pch = pch, asp = 1, ...)), warning = function(w) {message(w)})
        if (!is.null (k))
          graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
      }
      else if (type [1] == "som")
      {
        d = data.frame (Data = d)
        som = SOM (d, ...)
        if (is.null (k))
          graphics::plot (som, type = "mapping", labels = labels)
        else
        {
          graphics::plot (som, type = "mapping", col = col, labels = labels)
          graphics::legend (x = legendpos, legend = levels (k), fill = lcol, bty = "n")
        }
      }
    }
  }

#' Singular Value Decomposition
#'
#' Return the SVD decomposition.
#' @name SVD
#' @param x A numeric dataset (data.frame or matrix).
#' @param ndim The number of dimensions.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[base]{svd}}
#' @examples
#' require (datasets)
#' data (iris)
#' SVD (iris [, -5])
SVD <-
  function (x, ndim = min (nrow (x), ncol (x)), ...)
  {
    res = svd (x, nu = ndim, nv = ndim)
    ind = res$u %*% diag (res$d [1:ndim])
    rownames (ind) = rownames (x)
    colnames (ind) = paste ("Dim.", 1:ncol (ind))
    var = res$v %*% diag (res$d [1:ndim])
    rownames (var) = colnames (x)
    colnames (var) = paste ("Dim.", 1:ncol (var))
    proj = list (ind = ind, var = var)
    res = c (res, proj = list (proj))
    return (res)
  }

#' t-distributed Stochastic Neighbor Embedding
#'
#' Return the t-SNE dimensionality reduction.
#' @name TSNE
#' @param x A numeric dataset (data.frame or matrix).
#' @param perplexity Specification of the perplexity.
#' @param nstart How many random sets should be chosen?
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' @examples
#' require (datasets)
#' data (iris)
#' TSNE (iris [, -5])
TSNE <-
  function (x, perplexity = 30, nstart = 10, ...)
  {
    d = unique (x)
    res = Rtsne::Rtsne (d, perplexity = perplexity, ...)
    eval = utils::tail (res$itercosts, 1)
    for (i in 1:(nstart - 1))
    {
      tmp = Rtsne::Rtsne (d, perplexity = perplexity, ...)
      if (utils::tail (res$itercosts, 1) < eval)
      {
        res = tmp
        eval = utils::tail (res$itercosts, 1)
      }
    }
    rownames (res$Y) = rownames (d)
    colnames (res$Y) = paste ("Dim.", 1:ncol (res$Y))
    return (res)
  }
