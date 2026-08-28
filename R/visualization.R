#' @keywords internal
histogram <- function(d, main)
{
  histogram = graphics::hist (d, xlab = "", main = main, col = "#EB9292", ylab = "")
  graphics::mtext ("Frequency", side = 2, line = 3)
  xlim = c (min (histogram$breaks), max (histogram$breaks))
  dens = stats::density (d)
  # Scaled to the density's own maximum: a concentrated variable pushes it well above 1
  # (a N (0, 0.1) peaks near 4).
  ylim = c (0, max (dens$y))
  opar = graphics::par (new = TRUE)
  on.exit (graphics::par (opar))
  graphics::plot (dens, type = "l", col = "blue", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", xlim = xlim, ylim = ylim)
  graphics::axis (4)
  graphics::mtext ("Density", side = 4, line = 3)
}

#' Non-negative Matrix Factorization
#'
#' Return the NMF decomposition.
#' @name NMF
#' @param x A numeric dataset (data.frame or matrix).
#' @param rank Specification of the factorization rank.
#' @param nstart How many random sets should be chosen?
#' @param seed A specified seed for random number generation.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[NMF]{nmf}}
#' @examples
#' \dontrun{
#' install.packages ("BiocManager")
#' BiocManager::install ("Biobase")
#' install.packages ("NMF")
#' require (datasets)
#' data (iris)
#' NMF (iris [, -5])
#' }
NMF <-
  function (x, rank = 2, nstart = 10, seed = NULL, ...)
  {
    res = NULL
    if (requireNamespace ("NMF", quietly = TRUE))
    {
      setseed (seed)
      res = NMF::nmf (x, rank)
      eval = res@residuals
      for (i in seq_len (nstart - 1))
      {
        tmp = NMF::nmf (x, rank)
        if (tmp@residuals < eval)
        {
          res = tmp
          eval = res@residuals
        }
      }
      colnames (res@fit@W) = paste ("Dim.", 1:rank)
    }
    else
      message ("Package 'NMF' not installed!")
    return (res)
  }

#' @keywords internal
panel.hist <- function(x, ...)
{
  usr = graphics::par ("usr")
  on.exit (graphics::par (usr = usr))
  graphics::par (usr = c (usr [1:2], 0, 1.5))
  h = graphics::hist (x, plot = FALSE)
  breaks = h$breaks
  nB = length (breaks)
  y = h$counts
  y = y / max (y)
  graphics::rect (breaks [-nB], 0, breaks [-1], y, col = "grey")
}

#' @keywords internal
panel.blank <- function(x, y) {}

#' @keywords internal
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr = graphics::par ("usr")
  on.exit (graphics::par (usr = usr))
  graphics::par (usr = c (0, 1, 0, 1))
  r = stats::cor (x, y)
  txt = format (c (r, 0.123456789), digits  =digits) [1]
  txt = paste (prefix, txt, sep = "")
  if (missing (cex.cor))
    cex = 0.6 / graphics::strwidth (txt)
  graphics::text (0.5, 0.5, txt, cex = cex)
}

#' @keywords internal
# Shared low-level scatter plotting logic (color, point style, labels, legend, aspect ratio)
# used by every plotdata() branch that reduces the dataset to a 2-D projection
# (pairs/scatter, pca, cda, svd, nmf, tsne).
plotdata.scatter2d <-
  function (dd, d, k, col, pch, lcol, lpch, labels = FALSE, legendpos = "topleft", asp = 1, ...)
  {
    if (labels)
    {
      graphics::plot (dd, col = 0, asp = asp, ...)
      graphics::text (dd, row.names (d), col = col)
    }
    else
      graphics::plot (dd, col = col, pch = pch, asp = asp, ...)
    if (!is.null (k))
      graphics::legend (x = legendpos, legend = levels (k), pch = lpch, col = lcol, bty = "n")
  }

#' @keywords internal
# Splits the arguments of '...' between the method that computes a projection and the plot that
# draws it. plotdata() has to give some to one and some to the other, and gave all of them to
# both: 'perplexity', meant for TSNE(), reached graphics::plot() as well, which answered
# "\"perplexity\" is not a graphical parameter" once per element it drew. An argument goes to
# the method when it names one of its formals, and to the plot otherwise.
plotdata.route <-
  function (dots, fun)
  {
    named = names (dots)
    if (is.null (named))
      named = rep ("", length (dots))
    formal = setdiff (names (formals (fun)), "...")
    take = named %in% formal
    return (list (method = dots [take], plot = dots [!take]))
  }

#' @keywords internal
# How strongly each predictor relates to the target.
#
# A numeric target gives Pearson's correlation, sign included. A categorical one gives the
# correlation ratio -- the square root of the between-class share of the variance, which
# fseval.inertiaratio() computes -- because a correlation is not defined against more than two
# groups. On exactly two classes that ratio *is* |r| with the classes coded 0/1, so the sign of
# the gap between the two class means is put back: it says which way the variable moves.
plotdata.correlations <-
  function (d, target)
  {
    d = as.data.frame (d)
    numeric = sapply (d, is.numeric)
    if ((!any (numeric)) || is.null (target))
      return (NULL)
    d = d [, numeric, drop = FALSE]
    if (is.numeric (target))
      res = sapply (d, function (v) suppressWarnings (stats::cor (v, target)))
    else
    {
      target = factor (target)
      if (nlevels (target) < 2)
        return (NULL)
      res = sqrt (fseval.inertiaratio (d, target, vtype = "univariate"))
      if (nlevels (target) == 2)
        res = res * sign (sapply (d, function (v) diff (tapply (v, target, mean))))
    }
    # A constant variable correlates with nothing; cor() answers NA and warns.
    res [!is.finite (res)] = 0
    names (res) = colnames (d)
    return (res)
  }

#' @keywords internal
plotdata.correlationplot <-
  function (d, target, legendpos = "bottomright", ...)
  {
    res = plotdata.correlations (d, target)
    if (is.null (res) || (length (res) == 0))
    {
      message ("plotdata: type = \"correlation\" needs numeric variables and a target to ",
               "relate them to -- give one as the second argument, or leave a single ",
               "qualitative column in the data for it to be taken from.")
      return (invisible (NULL))
    }
    signed = any (res < 0)
    # Sorted by strength, and horizontally: barplot() draws the first bar at the bottom, so the
    # strongest has to come last for the picture to read from the top down.
    res = res [order (abs (res))]
    extra = max (graphics::strwidth (names (res), units = "figure") * 30)
    opar = graphics::par (mar = graphics::par ("mar") + c (0, extra, 0, 0))
    on.exit (graphics::par (opar))
    graphics::barplot (res, horiz = TRUE, las = 1, border = NA,
                       col = ifelse (res < 0, 2, 4),
                       xlim = c (if (signed) -1 else 0, 1),
                       xlab = if (is.numeric (target)) "Pearson correlation with the target"
                              else "correlation ratio (eta) with the target", ...)
    graphics::abline (v = 0, col = "grey40")
    if (signed)
      graphics::legend (legendpos, c ("positive", "negative"), fill = c (4, 2), bty = "n")
    return (invisible (res))
  }

#' @keywords internal
plotdata.pairsplot <-
  function (d, col, pch, asp, ...)
  {
    # pairs() draws one panel per pair of variables. Past a dozen or so they no longer fit, and
    # R stops on "figure margins too large" -- which says nothing about what to do instead.
    tryCatch (graphics::pairs (d, upper.panel = panel.cor, diag.panel = panel.hist,
                               pch = pch, col = col, asp = asp, ...),
              error = function (e)
                stop ("plotdata: a matrix of scatter plots needs one panel per pair of ",
                      "variables, and ", ncol (d), " of them do not fit on this device (",
                      conditionMessage (e), "). Enlarge the graphics window, keep a subset of ",
                      "the variables, or use a projection -- type = \"pca\", \"svd\" or ",
                      "\"tsne\" -- which shows them all on two axes.", call. = FALSE))
  }

#' @keywords internal
plotdata.boxplot.multi <-
  function (d, k, col, lcol, legendpos)
  {
    mini = min (d)
    maxi = max (d)
    names = colnames (d)
    v = 1
    at = NULL
    if (!is.null (k))
    {
      nbclusters = length (unique (col))
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

#' @keywords internal
plotdata.som.multi <-
  function (d, k, col, lcol, legendpos, labels = FALSE, ...)
  {
    d = data.frame (Data = d)
    som = do.call (SOM, c (list (d), plotdata.route (list (...), SOM)$method))
    if (is.null (k))
      graphics::plot (som, type = "mapping", labels = labels)
    else
    {
      graphics::plot (som, type = "mapping", col = col, labels = labels)
      graphics::legend (x = legendpos, legend = levels (k), fill = lcol, bty = "n")
    }
  }

#' @keywords internal
plotdata.heatmap.multi <-
  function (d, clustered = FALSE)
  {
    d = as.matrix (d)
    if (clustered)
      stats::heatmap (d, hclustfun = HCA, cexRow = 0.2 + 1 / log10 (nrow (d) * 10), cexCol = 0.2 + 1 / log10 (ncol (d) * 10))
    else
      stats::heatmap (d, Rowv = NA, Colv = NA, cexRow = 0.2 + 1 / log10 (nrow (d) * 10), cexCol = 0.2 + 1 / log10 (ncol (d) * 10))
  }

#' @keywords internal
plotdata.parallel.multi <-
  function (d, k)
  {
    col = NULL
    if (is.null (k))
    {
      n = nrow (d)
      col = grDevices::rainbow (n) [sample (n, n)]
    }
    else
      col = as.numeric (k) + 1
    MASS::parcoord (d, col = col)
  }

#' @keywords internal
# Shared multi-panel layout used by the "histogram", "barplot" and "pie" types:
# lays out one panel per column of d and calls FUN (column, name) for each.
plotdata.panels.multi <-
  function (d, FUN)
  {
    # 'nrow'/'ncol' as local names would shadow base::nrow and base::ncol for the rest of the
    # body -- harmless here only because ncol (d) is evaluated first.
    n = ncol (d)
    rows = round (sqrt (n))
    cols = ceiling (n / rows)
    graphics::layout (matrix (1:(rows * cols), ncol = cols, byrow = TRUE))
    on.exit (graphics::layout (1))
    for (i in 1:n)
      FUN (d [, i], colnames (d) [i])
  }

#' @keywords internal
# Handles plotdata() for a single vector (is.vector (d) == TRUE).
plotdata.vector <-
  function (d, k, type, legendpos, col, pch, lcol, lpch, ...)
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
      if (!is.null (k))
      {
        graphics::boxplot (d~k, ylim = c (mini, maxi), ylab = "", col = lcol, xaxt='n', xlab = "")
        graphics::legend (x = legendpos, legend = levels (k), fill = lcol, bty = "n")
      }
      else
        graphics::boxplot (d, ylim = c (mini, maxi), ylab = "", col = "grey", xaxt='n', xlab = "")
    }
    else if (type [1] == "histogram")
      histogram (d, "")
    else if (type [1] == "barplot")
      graphics::barplot (table (d), main = "", border = 0)
    else if (type [1] == "pie")
      graphics::pie (table (d), main = "", col = grDevices::colorRampPalette (c ("#E0E0FF", "#4F4FFF")) (nlevels (d)))
    else if (type [1] == "words")
      plotcloud (d, k = k, ...)
    else
      message ("Unavailable plot")
  }

#' @keywords internal
# Handles plotdata() for a dataset with more than one variable (matrix or data.frame).
plotdata.matrix <-
  function (d, k, type, legendpos, alpha, asp, labels, col, pch, lcol, lpch, tsne, nmf,
            target = NULL, ...)
  {
    if (type [1] == "correlation")
      plotdata.correlationplot (d, target, legendpos, ...)
    else if ((type [1] == "pairs") & (ncol (d) > 2))
      plotdata.pairsplot (d, col, pch, asp, ...)
    else if ((type [1] == "scatter") | (type [1] == "pairs"))
    {
      dd = d
      if (ncol (d) != 2)
        # scale.unit = FALSE is intentional here: plotdata() shows the data on their original
        # scale, unlike PCA() (whose scale.unit defaults to TRUE) which is meant for a proper
        # factorial analysis. Not exposed as a plotdata() argument -- see @param type above.
        dd = FactoMineR::PCA (d, scale.unit = FALSE, ncp = 2, graph = FALSE)$ind$coord [, 1:2]
      plotdata.scatter2d (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = asp, ...)
    }
    else if (type [1] == "boxplot")
      plotdata.boxplot.multi (d, k, col, lcol, legendpos)
    else if (type [1] == "pca")
    {
      # Same intentional scale.unit = FALSE as above.
      dd = FactoMineR::PCA (d, scale.unit = FALSE, ncp = 2, graph = FALSE)$ind$coord [, 1:2]
      plotdata.scatter2d (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = asp, ...)
    }
    else if (type [1] == "cda")
    {
      if (is.null (k))
        message ("Unavailable plot")
      else
      {
        proj = CDA (d, k)$proj
        if (ncol (proj) > 1)
          plotdata.scatter2d (proj [, 1:2], d, k, col, pch, lcol, lpch, labels, legendpos,
                              asp = asp, ...)
        else
        {
          # Two classes give a single canonical axis, and proj [, 1:2] then asked for a column
          # that does not exist. There is nothing to plot the axis against but the observation
          # index -- which is what plot.cda() does in the same situation, and what plotdata()
          # does for a dataset reduced to one variable. No aspect ratio: the two axes have no
          # common unit.
          dd = cbind (Index = seq_len (nrow (proj)), proj [, 1])
          colnames (dd) [2] = colnames (proj) [1]
          rownames (dd) = rownames (proj)
          plotdata.scatter2d (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = NA, ...)
        }
      }
    }
    else if (type [1] == "svd")
    {
      dd = SVD (d)$proj$ind [, 1:2]
      plotdata.scatter2d (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = asp, ...)
    }
    else if (type [1] == "nmf")
    {
      res = nmf
      if (is.null (res))
        res = NMF (d)
      dd = res@fit@W [, 1:2]
      plotdata.scatter2d (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = asp, ...)
    }
    else if (type [1] == "tsne")
    {
      res = tsne
      args = plotdata.route (list (...), TSNE)
      if (is.null (res))
        res = do.call (TSNE, c (list (d), args$method))
      dd = res$Y
      do.call (plotdata.scatter2d,
               c (list (dd, d, k, col, pch, lcol, lpch, labels, legendpos, asp = asp),
                  args$plot))
    }
    else if (type [1] == "som")
      plotdata.som.multi (d, k, col, lcol, legendpos, labels, ...)
    else if (type [1] == "heatmap")
      plotdata.heatmap.multi (d, clustered = FALSE)
    else if (type [1] == "heatmapc")
      plotdata.heatmap.multi (d, clustered = TRUE)
    else if (type [1] == "parallel")
      plotdata.parallel.multi (d, k)
    else if (type [1] == "histogram")
      plotdata.panels.multi (d, function (col, name) histogram (col, name))
    else if (type [1] == "barplot")
      plotdata.panels.multi (d, function (col, name) graphics::barplot (table (col), main = name, border = 0))
    else if (type [1] == "pie")
      plotdata.panels.multi (d, function (col, name) graphics::pie (table (col), main = name, col = grDevices::colorRampPalette (c ("#E0E0FF", "#4F4FFF")) (nlevels (col))))
    else
      # plotdata.vector() has always ended on such a message; without one here, an unsupported
      # type -- or a plain typo -- on a multi-variable dataset drew nothing and said nothing.
      # "words" lands here too: word clouds need a corpus of texts, see plotcloud().
      message (paste (type [1], ": unavailable plot for a dataset with several variables"))
  }

#' Advanced plot function
#'
#' Plot a dataset.
#'
#' \code{type = "correlation"} draws how strongly each variable relates to the target, sorted,
#' strongest at the top. \strong{Two different quantities, depending on the target.} Against a
#' numeric target it is Pearson's correlation \eqn{r}, sign included, and the axis says so.
#' Against a categorical one it is the \strong{correlation ratio} \eqn{\eta}, the square root
#' of the between-class share of the variance -- \emph{not} a Pearson coefficient computed on
#' class numbers, which would depend on the order the classes happen to be in and would mean
#' nothing beyond two classes. \eqn{\eta} lies in [0, 1], is defined for any number of classes
#' and does not depend on their coding. On exactly two classes \eqn{\eta} is \eqn{|r|} with the
#' classes coded 0/1, so the sign comes back and says which class the variable is larger in.
#'
#' The projections (\code{type = "pca"}, and the default \code{"scatter"}/\code{"pairs"} on
#' more than two variables) are computed on the data as they are, unscaled -- \code{plotdata}
#' shows a dataset, whereas \code{\link{PCA}} performs a factorial analysis and centres and
#' scales by default. So \code{plotdata (d, type = "pca")} and \code{plot (PCA (d))} differ
#' visibly when the variables have very different scales, and \code{\link{PCA}} is the one to
#' use for a properly scaled projection.
#' @name plotdata
#' @param d A numeric dataset (data.frame or matrix).
#' @param k The variable the observations are told apart by: they are coloured, grouped or,
#' for \code{type = "cda"}, discriminated by it. It is categorical, and cluster numbers do just
#' as well as a \code{factor}. Left \code{NULL}, it is taken from \code{d} when exactly one of
#' its columns is qualitative, or from \code{target} when that one is categorical.
#' @param target The variable to be explained, read by \code{type = "correlation"} only. Unlike
#' \code{k} it may be continuous. Give one or the other, the way \code{\link[stats]{cutree}}
#' takes either \code{k} or \code{h}: a categorical \code{target} also serves as \code{k},
#' and \code{k} serves as \code{target} when none is given. A continuous \code{target} leaves
#' the observations uncoloured, having no groups to offer.
#' @param type The type of graphic to be plotted. See the Details section on the projections.
#' @param legendpos Position of the legend
#' @param alpha Opacity of the plotted points, from 0 (invisible) to 255 (opaque). Useful on
#' dense scatter plots, where points would otherwise hide each other. The legend stays
#' opaque.
#' @param asp Aspect ratio: 1 (the default) makes one unit as long on both axes, \code{NA}
#' lets them scale independently.
#' @param labels Indicates whether or not labels (row names) should be showned on the (scatter) plot.
#' @param tsne A precomputed \code{\link{TSNE}} result. When \code{type = "tsne"}, providing this avoids
#' recomputing the (randomized, potentially costly) t-SNE embedding on every call; if \code{NULL} (default),
#' it is computed internally as before.
#' @param nmf A precomputed \code{\link{NMF}} result. When \code{type = "nmf"}, providing this avoids
#' recomputing the (randomized, potentially costly) NMF decomposition on every call; if \code{NULL} (default),
#' it is computed internally as before.
#' @param ... Other parameters.
#' @export
#' @examples
#' require (datasets)
#' data (iris)
#' # Without classification
#' plotdata (iris [, -5]) # Default (pairs)
#' # With classification
#' plotdata (iris [, -5], iris [, 5]) # Default (pairs)
#' plotdata (iris, 5) # Column number
#' plotdata (iris) # Automatic detection of the classification (if only one factor column)
#' plotdata (iris, type = "scatter") # Scatter plot (PCA axis)
#' plotdata (iris, type = "parallel") # Parallel coordinates
#' plotdata (iris, type = "boxplot") # Boxplot
#' plotdata (iris, type = "histogram") # Histograms
#' plotdata (iris, type = "heatmap") # Heatmap
#' plotdata (iris, type = "heatmapc") # Heatmap (and hierarchalcal clustering)
#' plotdata (iris, type = "pca") # Scatter plot (PCA axis)
#' plotdata (iris, type = "cda") # Scatter plot (CDA axis)
#' plotdata (iris, type = "svd") # Scatter plot (SVD axis)
#' plotdata (iris, type = "som") # Kohonen map
#' # With only one variable
#' plotdata (iris [, 1], iris [, 5]) # Default (data vs. index)
#' plotdata (iris [, 1], iris [, 5], type = "scatter") # Scatter plot (data vs. index)
#' plotdata (iris [, 1], iris [, 5], type = "boxplot") # Boxplot
#' # With two variables
#' plotdata (iris [, 3:4], iris [, 5]) # Default (scatter plot)
#' plotdata (iris [, 3:4], iris [, 5], type = "scatter") # Scatter plot
#' data (titanic)
#' plotdata (titanic, type = "barplot") # Barplots
#' plotdata (titanic, type = "pie") # Pie charts
#' \dontrun{
#' # Reusing a previously computed t-SNE embedding instead of recomputing it
#' res = TSNE (iris [, -5])
#' plotdata (iris [, -5], iris [, 5], type = "tsne", tsne = res)
#' }
plotdata <-
  function (d, k = NULL, target = NULL,
            type = c ("pairs", "scatter", "parallel", "boxplot", "histogram", "barplot", "pie", "heatmap", "heatmapc", "correlation", "pca", "cda", "svd", "nmf", "tsne", "som", "words"),
            legendpos = "topleft", alpha = 200, asp = 1, labels = FALSE, tsne = NULL, nmf = NULL, ...)
  {
    # The bars of a correlation plot are sorted by strength, so the two bottom corners are the
    # ones with room; "topleft", which suits a scatter plot, sits on the longest bar. Only the
    # default is overridden -- an explicit legendpos is honoured.
    if (missing (legendpos) && (type [1] == "correlation"))
      legendpos = "bottomright"
    factors = NULL
    if (is.factor (d))
      factors = TRUE
    else if (is.vector (d))
      factors = FALSE
    else
      factors = sapply (as.data.frame (d), is.factor)
    if ((type [1] == "barplot") || (type [1] == "pie"))
      d = d [, factors]
    else
    {
      if ((is.null (k)) && (sum (factors) == 1))
        k = d [, factors]
      else if ((length (k) == 1) && (factors [k]))
        k = d [, k]
      if (sum (factors) > 0)
        d = d [, !factors]
    }
    # 'k' and 'target' fill in for each other, the way cutree() takes either 'k' or 'h': 'k'
    # tells the observations apart (colour, groups, discrimination) and has to be categorical;
    # 'target' is the variable to be explained by type = "correlation" and may be continuous.
    # A categorical target serves as 'k' as well; a continuous one cannot, and leaves the
    # observations uncoloured.
    if (is.null (target))
      target = k
    else if (is.null (k) && (!is.numeric (target)))
      k = target
    col = 1
    pch = 1
    # 'k' becomes a factor only where it is read as a class -- to colour, to group, or to
    # discriminate on. type = "correlation" reads it as the target instead, on its own scale:
    # coercing a numeric one would build as many levels as there are observations, and as many
    # colours, all of them thrown away.
    if ((!is.null (k)) && (type [1] != "correlation"))
    {
      if (!is.factor (k))
        k = factor (k)
      pch = as.numeric (k) + 1
      col = pch
    }
    # The legend keeps the opaque colours; only the plotted points get the alpha channel.
    lcol = sort (unique (col))
    lpch = sort (unique (pch))
    # 'alpha' was documented ("Color opacity (0-255)") and passed down to plotdata.matrix(),
    # but never used anywhere -- as was addalpha(), the helper written for exactly this.
    # Applying it makes dense scatter plots readable, which is what the parameter is for;
    # alpha = 255 restores fully opaque points.
    if ((!is.null (alpha)) && (alpha < 255))
      col = addalpha (col, alpha)
    if (length (d) == 0)
      message ("Unavailable plot")
    else if (is.vector (d))
      plotdata.vector (d, k, type, legendpos, col, pch, lcol, lpch, ...)
    else
      plotdata.matrix (d, k, type, legendpos, alpha, asp, labels, col, pch, lcol, lpch, tsne,
                       nmf, target, ...)
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
    maxdim = min (nrow (x), ncol (x))
    if ((!is.numeric (ndim)) || (length (ndim) != 1) || is.na (ndim) ||
        (ndim < 1) || (ndim > maxdim) || (ndim != round (ndim)))
      stop ("SVD: 'ndim' must be a single whole number between 1 and ", maxdim,
            " (the smaller dimension of 'x'). Got: ", ndim, ".")
    res = svd (x, nu = ndim, nv = ndim)
    # 'nrow = ndim' is not optional: diag() of a *single* number builds an identity matrix of
    # that size instead of a 1 x 1 matrix holding it.
    d = diag (res$d [1:ndim], nrow = ndim)
    ind = res$u %*% d
    rownames (ind) = rownames (x)
    colnames (ind) = paste ("Dim.", 1:ncol (ind))
    var = res$v %*% d
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
#' @param nstart How many random sets should be chosen? The embedding with the lowest final
#' cost is kept.
#' @param seed A specified seed for random number generation.
#' @param ... Other parameters.
#' @return The \code{\link[Rtsne]{Rtsne}} result. \code{Rtsne} requires distinct observations,
#' so duplicated rows of \code{x} are removed before fitting; the returned \code{Y} (and
#' \code{costs}) are then expanded back to \code{nrow (x)} rows, in the order of \code{x}, so
#' that they can be used directly alongside the original dataset -- duplicated observations
#' simply share the same coordinates.
#' @export
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' @examples
#' require (datasets)
#' data (iris)
#' TSNE (iris [, -5])
TSNE <-
  function (x, perplexity = 30, nstart = 10, seed = NULL, ...)
  {
    # Rtsne rejects duplicated rows, so the embedding is fitted on the distinct ones and
    # expanded back afterwards: res$Y must have one row per row of x, in the same order, or
    # callers such as plotdata (d, k, type = "tsne") mismatch points and colours.
    keys = do.call (paste, c (as.data.frame (x), sep = "\r"))
    keep = !duplicated (keys)
    d = x [keep, , drop = FALSE]
    setseed (seed)
    res = Rtsne::Rtsne (d, perplexity = perplexity, ...)
    eval = utils::tail (res$itercosts, 1)
    for (i in seq_len (nstart - 1))
    {
      tmp = Rtsne::Rtsne (d, perplexity = perplexity, ...)
      # 'tmp', not 'res': comparing the current best with its own cost is always FALSE, so
      # the extra nstart - 1 embeddings were computed and then thrown away.
      if (utils::tail (tmp$itercosts, 1) < eval)
      {
        res = tmp
        eval = utils::tail (res$itercosts, 1)
      }
    }
    back = match (keys, keys [keep])
    res$Y = res$Y [back, , drop = FALSE]
    if (!is.null (res$costs))
      res$costs = res$costs [back]
    rownames (res$Y) = rownames (as.data.frame (x))
    colnames (res$Y) = paste ("Dim.", 1:ncol (res$Y))
    return (res)
  }

