#' DBSCAN model
#'
#' This class contains the model obtained by the DBSCAN method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{cluster}}{A vector of integers indicating the cluster to which each point is allocated.}
#'   \item{\code{eps}}{Reachability distance (parameter).}
#'   \item{\code{MinPts}}{Reachability minimum no. of points (parameter).}
#'   \item{\code{isseed}}{A logical vector indicating whether a point is a seed (not border, not noise).}
#'   \item{\code{data}}{The dataset that has been used to fit the map (as a \code{matrix}).}
#' }
#' @name dbs-class
#' @seealso \code{\link{DBSCAN}}
NULL

#' Expectation-Maximization model
#'
#' This class contains the model obtained by the EM method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{modelName}}{A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.}
#'   \item{\code{prior}}{Specification of a conjugate prior on the means and variances.}
#'   \item{\code{n}}{The number of observations in the dataset.}
#'   \item{\code{d}}{The number of variables in the dataset.}
#'   \item{\code{G}}{The number of components of the mixture.}
#'   \item{\code{z}}{A matrix whose \code{[i,k]}th entry is the conditional probability of the ith observation belonging to the kth component of the mixture.}
#'   \item{\code{parameters}}{A names list giving the parameters of the model.}
#'   \item{\code{control}}{A list of control parameters for EM.}
#'   \item{\code{loglik}}{The log likelihood for the data in the mixture model.}
#'   \item{\code{cluster}}{A vector of integers (from \code{1:k}) indicating the cluster to which each point is allocated.}
#' }
#' @name em-class
#' @seealso \code{\link{EM}}, \code{\link[mclust]{mclustModelNames}}
NULL

#' MeanShift model
#'
#' This class contains the model obtained by the MEANSHIFT method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{cluster}}{A vector of integers indicating the cluster to which each point is allocated.}
#'   \item{\code{value}}{A vector or matrix containing the location of the classified local maxima in the support.}
#'   \item{\code{data}}{The leaning set.}
#'   \item{\code{kernel}}{A string indicating the kernel associated with the kernel density estimate that the mean shift is optimizing over.}
#'   \item{\code{bandwidth}}{Used in the kernel density estimate for steepest ascent classification.}
#'   \item{\code{alpha}}{A scalar tuning parameter for normal kernels.}
#'   \item{\code{iterations}}{The number of iterations to perform mean shift.}
#'   \item{\code{epsilon}}{A scalar used to determine when to terminate the iteration of an individual query point.}
#'   \item{\code{epsilonCluster}}{A scalar used to determine the minimum distance between distinct clusters.}
#' }
#' @name meanshift-class
#' @seealso \code{\link{MEANSHIFT}}
NULL

#' Self-Organizing Maps model
#'
#' This class contains the model obtained by the SOM method.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{som}}{An object of class \code{kohonen} representing the fitted map.}
#'   \item{\code{nodes}}{A \code{vector} of integer indicating the cluster to which each node is allocated.}
#'   \item{\code{cluster}}{A \code{vector} of integer indicating the cluster to which each observation is allocated.}
#'   \item{\code{data}}{The dataset that has been used to fit the map (as a \code{matrix}).}
#' }
#' @name som-class
#' @seealso \code{\link{plot.som}}, \code{\link{SOM}}, \code{\link[kohonen]{som}}
NULL

#' Spectral clustering model
#'
#' This class contains the model obtained by Spectral clustering.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{cluster}}{A \code{vector} of integer indicating the cluster to which each observation is allocated.}
#'   \item{\code{proj}}{The projection of the dataset in the spectral space.}
#'   \item{\code{centers}}{The cluster centers (on the spectral space).}
#'   \item{\code{data}}{The dataset that has been used to fit the model (as a \code{matrix}).}
#' }
#' @name spectral-class
#' @seealso \code{\link{SPECTRAL}}, \code{\link{predict.spectral}}
NULL

#' @keywords internal
pca.project2d <-
  function (d)
  {
    pca = FactoMineR::PCA (d, scale.unit = FALSE, graph = FALSE, ncp = 2)
    coord = pca$ind$coord
    row.names (coord) = row.names (d)
    # Projects any point of the original space onto exactly the same two axes as 'coord'
    # (FactoMineR's own transform: centre, divide by ecart.type -- all ones here, since
    # scale.unit = FALSE -- then rotate by svd$V).
    #
    # Shared by scatterplot() and plot.som(), which both have to place cluster centres on the
    # very axes the point cloud was projected onto.
    project =
      function (m)
        sweep (sweep (as.matrix (m), 2, pca$call$centre, "-"), 2, pca$call$ecart.type, "/") %*% pca$svd$V
    xlab = paste ("Prin. 1 (", round (pca$eig [1, 2], 2), " %)", sep = "")
    ylab = paste ("Prin. 2 (", round (pca$eig [2, 2], 2), " %)", sep = "")
    return (list (pca = pca, coord = coord, project = project, xlab = xlab, ylab = ylab))
  }

#' @keywords internal
# Agreement between two binary indicator vectors. Goes through align.labels()
# (classification.R) for the same reason as eval.accuracy(): table() drops the values it
# does not observe, so a cluster indicator that happens to be constant (all TRUE or all
# FALSE) produced a 1 x 2 table whose diagonal read the wrong cells.
accuracy0 <-
  function (clus, gt)
  {
    a = align.labels (clus, gt)
    return (sum (a$predictions == a$gt, na.rm = TRUE) / length (a$gt))
  }

#' @keywords internal
# Best agreement between one cluster indicator and any single class of 'gt'.
#
# Iterates over the values 'gt' actually takes -- they need be neither contiguous nor start at
# 1 -- and drops the comparisons that come out NaN.
bestoverclasses <-
  function (FUN, clus, gt)
  {
    res = sapply (sort (unique (gt)), function (i) FUN (clus, gt == i))
    if (all (is.na (res)))
      return (NA)
    return (max (res, na.rm = TRUE))
  }

#' @keywords internal
accuracy1 <-
  function (clus, gt)
  {
    if (sum (clus) == 0)
      return (NA)
    return (bestoverclasses (accuracy0, clus, gt))
  }

#' Clustering Box Plots
#'
#' Produce a box-and-whisker plot for clustering results.
#' @name boxclus
#' @inheritParams scatterplot
#' @param legendpos Position of the legend
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[graphics]{boxplot}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' boxclus (iris [, -5], km$cluster)
boxclus <-
  function (d, clusters, legendpos = "topleft", ...)
  {
    # 'clusters' is documented as a vector *or a factor*, but the body treated it as numeric
    # throughout -- min(), 1 + clusters -- so a factor (a ground truth, the classes of a
    # dataset) failed on "'min' not meaningful for factors" where scatterplot() accepted it.
    # Both go through the same helper now.
    coded = cluster.codes (clusters)
    clusters = coded$codes
    mini = min (d)
    maxi = max (d)
    nbclusters = length (unique (clusters))
    names = colnames (d)
    at = (0:(ncol (d) - 1)) * nbclusters + 1 + (nbclusters - 1) / 2
    v = (1:(ncol (d) - 1)) * nbclusters + .5
    d = utils::stack (d)
    d$cluster = clusters
    # The boxes are coloured 1 + code, as everywhere else in the package, so that the legend
    # below matches them and noise (cluster 0) comes out black.
    legend = cluster.legend (clusters, coded$names)
    graphics::boxplot (values~cluster+ind, d, ylim = c (mini, maxi), ylab = "",
                       col = legend$col, xaxt = 'n', xlab = "")
    graphics::axis (side = 1, at = at, labels = names, lwd.ticks = FALSE, lwd = 0)
    graphics::abline (v = v, lty = 2, col = "grey")
    graphics::legend (x = legendpos, legend = legend$labels, fill = legend$col, bty = "n")
  }

#' @keywords internal
# Lookup table for the compare.* dispatch used by compare(). Built lazily (as a function
# rather than a static list) so it does not depend on file collation order.
compare.functions <-
  function ()
  {
    list (accuracy = compare.accuracy,
          jaccard = compare.jaccard,
          kappa = compare.kappa)
  }

#' Comparison of two sets of clusters
#'
#' Comparison of two sets of clusters
#' @name compare
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @param eval The evaluation criterion.
#' @param comp How the two partitions are compared: \code{"max"} (each cluster matched with the
#' class it agrees with most, averaged over clusters), \code{"cluster"} (the same scores, not
#' averaged) or \code{"pairwise"} (every pair of observations, labels ignored). In
#' \code{"pairwise"} mode the three criteria are three classical indices built on the same pair
#' counts: \code{"accuracy"} is the Rand index, \code{"jaccard"} the Jaccard index on pairs,
#' and \code{"kappa"} Cohen's kappa on the fourfold table of pair agreements (which, by
#' Warrens (2008), is also the adjusted Rand index). See \code{\link{compare.accuracy}},
#' \code{\link{compare.jaccard}} and \code{\link{compare.kappa}}.
#' @return A numeric value indicating how much the two sets of clusters are similar.
#' @export
#' @seealso \code{\link{compare.accuracy}}, \code{\link{compare.jaccard}}, \code{\link{compare.kappa}}, \code{\link{intern}}, \code{\link{stability}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' compare (km$cluster, iris [, 5])
#' \dontrun{
#' compare (km$cluster, iris [, 5], eval = c ("accuracy", "kappa"), comp = "pairwise")
#' }
compare <-
  function (clus, gt, eval = "accuracy", comp = c ("max", "pairwise", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    comp = match.arg (comp [1], c ("max", "pairwise", "cluster"))
    funs = compare.functions ()
    res = lapply (eval, function (e)
    {
      fun = funs [[e]]
      if (is.null (fun))
        stop ("compare: unknown evaluation criterion '", e, "'. Available criteria: ",
              paste (names (funs), collapse = ", "))
      return (fun (clus, gt, comp))
    })
    if (comp == "cluster")
    {
      # One value per cluster and per criterion. Naming the whole thing after the criteria --
      # which is what the other two modes need -- overwrote the cluster names and padded with
      # NA: two criteria on three clusters came back as six values named
      # "jaccard", "accuracy", NA, NA, NA, NA. A matrix of criteria by clusters, as intern()
      # returns in the same situation.
      out = do.call (rbind, res)
      rownames (out) = eval
      if (length (eval) == 1)
      {
        # The cluster names are restored by hand: on a single criterion *and* a single cluster,
        # 'out' is a 1x1 matrix and out [1, ] drops to a bare, unnamed scalar. stability()
        # matches its runs to the reference partition by those names, so it then found none and
        # averaged NA -- a degenerate clustering (mean shift returning one cluster) came back
        # as NaN instead of 1.
        one = out [1, ]
        names (one) = colnames (out)
        return (one)
      }
      return (out)
    }
    out = unlist (res)
    names (out) = eval
    return (out)
  }

#' Comparison of two sets of clusters, using accuracy
#'
#' Comparison of two sets of clusters, using accuracy
#' @name compare.accuracy
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @param comp How the two partitions are compared. \code{"max"} matches each extracted
#' cluster with the class it agrees with most and averages the per-cluster scores, weighted by
#' cluster size; \code{"cluster"} returns those scores instead of averaging them;
#' \code{"pairwise"} ignores the labels and looks at every pair of observations, asking whether
#' the two partitions agree on grouping it or separating it.
#' In \code{"pairwise"} mode this function is the \strong{Rand index}: the proportion of pairs
#' the two partitions agree on, counting both those they group and those they keep apart. See
#' \code{\link{compare.jaccard}} and \code{\link{compare.kappa}} for the two other readings
#' of the same pair counts.
#' @return A numeric value indicating how much the two sets of clusters are similar.
#' @export
#' @seealso \code{\link{compare.jaccard}}, \code{\link{compare.kappa}}, \code{\link{compare}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' compare.accuracy (km$cluster, iris [, 5])
compare.accuracy <-
  function (clus, gt, comp = c ("max", "pairwise", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    comp = match.arg (comp [1], c ("max", "pairwise", "cluster"))
    # as.numeric() on a character ground truth -- read.table() gives one whenever the class
    # column is not converted to a factor -- returns NA for every observation, and the pair
    # counts below then came out as a spurious 1. cluster.codes() maps categories to codes.
    kk1 = cluster.codes (clus)$codes
    kk2 = cluster.codes (gt)$codes
    res = 0
    if (comp == "pairwise")
    {
      res = pairwise.rand (kk1, kk2)
    }
    else
    {
      res = NULL
      # The cluster values actually present, not min:max -- otherwise a partition whose ids
      # have a gap (cluster 2 empty, say) produced one entry per *possible* id while
      # table (kk1) only counts the existing ones, and weighted.mean() below failed with
      # "'x' and 'w' must have the same length".
      clusters = sort (unique (kk1))
      for (i in clusters)
        res = c (res, accuracy1 (kk1 == i, kk2))
      if (comp == "max")
        res = stats::weighted.mean (res, table (kk1), na.rm = TRUE)
      else
        names (res) = paste ("Cluster", clusters)
    }
    return (res)
  }

#' Comparison of two sets of clusters, using Jaccard index
#'
#' Comparison of two sets of clusters, using Jaccard index
#' @name compare.jaccard
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @inheritParams compare.accuracy
#' @section The pairwise index: the \strong{Jaccard index} on pairs -- the pairs both
#' partitions group together, over the pairs at least one of them groups. Unlike the Rand index
#' of \code{\link{compare.accuracy}} it ignores the pairs both keep apart, a cell that
#' dominates as soon as there are many clusters: 150 singletons out of 150 observations score
#' 0.67 with the Rand index and 0 here.
#' @return A numeric value indicating how much the two sets of clusters are similar.
#' @export
#' @seealso \code{\link{compare.accuracy}}, \code{\link{compare.kappa}}, \code{\link{compare}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' compare.jaccard (km$cluster, iris [, 5])
compare.jaccard <-
  function (clus, gt, comp = c ("max", "pairwise", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    comp = match.arg (comp [1], c ("max", "pairwise", "cluster"))
    # as.numeric() on a character ground truth -- read.table() gives one whenever the class
    # column is not converted to a factor -- returns NA for every observation, and the pair
    # counts below then came out as a spurious 1. cluster.codes() maps categories to codes.
    kk1 = cluster.codes (clus)$codes
    kk2 = cluster.codes (gt)$codes
    res = 0
    if (comp == "pairwise")
    {
      res = pairwise.jaccard (kk1, kk2)
    }
    else
    {
      res = NULL
      # The cluster values actually present, not min:max -- otherwise a partition whose ids
      # have a gap (cluster 2 empty, say) produced one entry per *possible* id while
      # table (kk1) only counts the existing ones, and weighted.mean() below failed with
      # "'x' and 'w' must have the same length".
      clusters = sort (unique (kk1))
      for (i in clusters)
        res = c (res, jaccard1 (kk1 == i, kk2))
      if (comp == "max")
        res = stats::weighted.mean (res, table (kk1), na.rm = TRUE)
      else
        names (res) = paste ("Cluster", clusters)
    }
    return (res)
  }

#' Comparison of two sets of clusters, using kappa
#'
#' Comparison of two sets of clusters, using kappa
#' @name compare.kappa
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @inheritParams compare.accuracy
#' @section The pairwise index: \strong{Cohen's kappa} on the fourfold table of pair
#' agreements, i.e. the Rand index of \code{\link{compare.accuracy}} corrected for the
#' agreement expected by chance. That quantity is also the Hubert-Arabie \emph{adjusted Rand
#' index} -- a theorem of Warrens (2008), not a substitution.
#' @references Warrens, M.J. (2008). On the Equivalence of Cohen's Kappa and the Hubert-Arabie
#' Adjusted Rand Index. \emph{Journal of Classification}, 25(2), 177-183.
#' @return A numeric value indicating how much the two sets of clusters are similar.
#' @export
#' @seealso \code{\link{compare.accuracy}}, \code{\link{compare.jaccard}}, \code{\link{compare}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' compare.kappa (km$cluster, iris [, 5])
compare.kappa <-
  function (clus, gt, comp = c ("max", "pairwise", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    comp = match.arg (comp [1], c ("max", "pairwise", "cluster"))
    # as.numeric() on a character ground truth -- read.table() gives one whenever the class
    # column is not converted to a factor -- returns NA for every observation, and the pair
    # counts below then came out as a spurious 1. cluster.codes() maps categories to codes.
    kk1 = cluster.codes (clus)$codes
    kk2 = cluster.codes (gt)$codes
    res = 0
    if (comp == "pairwise")
    {
      res = pairwise.kappa (kk1, kk2)
    }
    else
    {
      res = NULL
      # The cluster values actually present, not min:max -- otherwise a partition whose ids
      # have a gap (cluster 2 empty, say) produced one entry per *possible* id while
      # table (kk1) only counts the existing ones, and weighted.mean() below failed with
      # "'x' and 'w' must have the same length".
      clusters = sort (unique (kk1))
      for (i in clusters)
        res = c (res, kappa1 (kk1 == i, kk2))
      if (comp == "max")
        res = stats::weighted.mean (res, table (kk1), na.rm = TRUE)
      else
        names (res) = paste ("Cluster", clusters)
    }
    return (res)
  }

#' @keywords internal
# How the two partitions classify each of the n (n - 1) / 2 pairs of observations: together
# in both, together in the first only, together in the second only, separated in both.
#
# They follow from the k1 x k2 contingency table in O(n + k1 k2), without ever materialising
# one row per pair -- which would be an n (n - 1) / 2 by n matrix, 29.8 GB for n = 2000.
pairwise.counts <-
  function (k1, k2)
  {
    tt = table (k1, k2)
    total = choose (length (k1), 2)
    both = sum (choose (tt, 2))
    only1 = sum (choose (rowSums (tt), 2)) - both
    only2 = sum (choose (colSums (tt), 2)) - both
    return (list (both = both, only1 = only1, only2 = only2,
                  neither = total - both - only1 - only2, total = total))
  }

#' @keywords internal
# Proportion of pairs the two partitions agree on (the Rand index).
pairwise.rand <-
  function (k1, k2)
  {
    p = pairwise.counts (k1, k2)
    return ((p$both + p$neither) / p$total)
  }

#' @keywords internal
# The Jaccard index between the two "is this pair grouped together?" indicators: the pairs the
# two partitions agree to group, over the pairs at least one of them groups. Unlike the Rand
# index above it ignores 'neither', i.e. the pairs both partitions keep apart -- which is what
# makes it informative when the number of clusters is large, since that cell then dominates.
pairwise.jaccard <-
  function (k1, k2)
  {
    p = pairwise.counts (k1, k2)
    denominator = p$both + p$only1 + p$only2
    if (denominator == 0)
      return (1)
    return (p$both / denominator)
  }

#' @keywords internal
# Cohen's kappa between the two "is this pair grouped together?" indicators -- what
# irr::kappa2 (cbind (p1, p2)) computed on the explicit pair vectors, and still verified
# against it by the test suite.
#
# By Warrens (2008), Journal of Classification 25(2):177-183, Cohen's kappa computed on this
# particular fourfold table is also the Hubert-Arabie adjusted Rand index. That is an identity
# between two definitions, not a change of statistic: what is computed here is, and remains,
# Cohen's kappa.
pairwise.kappa <-
  function (k1, k2)
  {
    p = pairwise.counts (k1, k2)
    po = (p$both + p$neither) / p$total
    pe = ((p$both + p$only1) * (p$both + p$only2) +
          (p$only2 + p$neither) * (p$only1 + p$neither)) / p$total^2
    return ((po - pe) / (1 - pe))
  }

#' @keywords internal
# A default value for DBSCAN's 'eps': the knee of the sorted minpts-distance curve, i.e. the
# point furthest from the chord joining its two ends -- the construction distplot() invites
# the user to do by eye. Both axes are rescaled to [0, 1] first, since they have no common
# unit.
dbscan.eps <-
  function (d, minpts)
  {
    if (minpts + 1 > nrow (as.matrix (d)))
      stop ("DBSCAN: 'minpts' (", minpts, ") is too large for a dataset of ",
            nrow (as.matrix (d)), " observations.")
    kdist = sort (kdistances (d, minpts), decreasing = TRUE)
    n = length (kdist)
    if ((n < 3) || (diff (range (kdist)) == 0))
      return (stats::median (kdist))
    return (unname (kdist [knee.index (kdist)]))
  }

#' @keywords internal
# The knee of a monotone curve: the point furthest from the chord joining its two ends -- the
# construction one does by eye on a k-distance plot or on an elbow plot. Both axes are
# rescaled to [0, 1] first, since they have no common unit. Returns the index of that point.
knee.index <-
  function (y)
  {
    n = length (y)
    if ((n < 3) || (diff (range (y)) == 0))
      return (max (1, round (n / 2)))
    x = (seq_len (n) - 1) / (n - 1)
    y = (y - min (y)) / diff (range (y))
    x1 = x [1]; y1 = y [1]; x2 = x [n]; y2 = y [n]
    gap = abs ((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / sqrt ((y2 - y1)^2 + (x2 - x1)^2)
    return (which.max (gap))
  }

#' DBSCAN clustering method
#'
#' Run the DBSCAN algorithm for clustering.
#' @name DBSCAN
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param minpts Reachability minimum no. of points.
#' @param eps Reachability distance. If \code{NULL} (the default), it is set to the knee of
#' the sorted \code{minpts}-distance curve -- the construction \code{\link{distplot}} invites
#' you to do by eye -- and the value used is reported in a message. This is only a starting
#' point: \code{eps} is the parameter DBSCAN is most sensitive to, and it is worth looking at
#' \code{distplot (minpts, d)} before settling on one.
#' @param graph A logical indicating whether or not a graphic should be plotted (the
#' \code{minpts}-distance curve used to choose \code{eps}, when \code{eps} is not given).
#' @param ... Other parameters.
#' @return A clustering model obtained by DBSCAN.
#' @export
#' @seealso \code{\link[fpc]{dbscan}}, \code{\link{dbs-class}}, \code{\link{distplot}}, \code{\link{predict.dbs}}
#' @examples
#' require (datasets)
#' data (iris)
#' DBSCAN (iris [, -5], minpts = 5, eps = 1)
DBSCAN <-
  function (d, minpts = 5, eps = NULL, graph = FALSE, ...)
  {
    if (is.null (eps))
    {
      eps = dbscan.eps (d, minpts)
      message ("DBSCAN: 'eps' not given, using ", signif (eps, 4),
               " (knee of the ", minpts, "-distance curve). Check it with distplot (",
               minpts, ", d) and pass 'eps' explicitly if another value suits better.")
      if (graph)
        distplot (minpts, d, h = eps)
    }
    res = fpc::dbscan (d, MinPts = minpts, eps = eps)
    res = c (res, list (data = d))
    class (res) = "dbs"
    return (res)
  }

#' Plot a k-distance graphic
#'
#' Plot the distance to the k's nearest neighbours of each object in decreasing order. Mostly used to determine the \code{eps} parameter for the \code{\link[fpc]{dbscan}} function.
#' @name distplot
#' @param k The \code{k} parameter.
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param h The y-coordinate at which a horizontal line should be drawn.
#' @export
#' @seealso \code{\link{DBSCAN}}, \code{\link[fpc]{dbscan}}
#' @examples
#' require (datasets)
#' data (iris)
#' distplot (5, iris [, -5], h = .65)
distplot <-
  function (k, d, h = -1)
  {
    Kdistance = sort (kdistances (d, k), decreasing = TRUE)
    graphics::plot (Kdistance, t = "l", xaxt = "n", xlab = "Objects", ylab = paste (k, "-distance"), main = "", xaxs = "i", yaxs = "i", col = "darkblue")
    if (h > 0)
    {
      graphics::abline (h = h, lty = 2, col = "blue")
      y = h
      x = 3 * nrow (d) / 4
      graphics::text (x, y, bquote (epsilon == .(h)), pos = 3, col = "blue")
    }
  }

#' Expectation-Maximization clustering method
#'
#' Run the EM algorithm for clustering.
#' @name EM
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k Either an integer (the number of clusters) or a (\code{vector}) indicating the cluster to which each point is initially allocated.
#' @param model A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @param seed A specified seed for random number generation (used only for the default k-means initialization).
#' @param ... Other parameters.
#' @return A clustering model obtained by EM.
#' @export
#' @seealso \code{\link[mclust]{em}}, \code{\link[mclust]{mstep}}, \code{\link[mclust]{mclustModelNames}}
#' @examples
#' require (datasets)
#' data (iris)
#' EM (iris [, -5], 3) # Default initialization
#' km = KMEANS (iris [, -5], k = 3)
#' EM (iris [, -5], km$cluster) # Initialization with another clustering method
EM <-
  function (d, k, model = "VVV", seed = NULL, ...)
  {
    clusters = k
    if (length (clusters) == 1)
    {
      setseed (seed)
      clusters = stats::kmeans (d, clusters, nstart = 10)$cluster
    }
    z = mclust::unmap (clusters)
    p = mclust::mstep (data = d, modelName = model, z = z)
    res = mclust::em (data = d, modelName = model, parameters = p$parameters)
    # mclust::em() has not returned 'loglik' in every release. One E-step on the fitted
    # parameters recovers it at negligible cost, so print.em() does not depend on the version
    # of mclust installed.
    if (is.null (res$loglik) || !is.finite (res$loglik))
      res$loglik = tryCatch (mclust::estep (data = d, modelName = model,
                                            parameters = res$parameters)$loglik,
                             error = function (e) NULL)
    res = c (res, list (cluster = apply (res$z, 1, which.max)))
    class (res) = "em"
    return (res)
  }

#' @keywords internal
# The linkage names agnes uses, mapped to the ones hclust uses. agnes's "ward" is hclust's
# "ward.D2": on every dataset tried the two give the same heights, to machine precision.
hca.linkage <-
  function (method)
  {
    map = c (ward = "ward.D2", single = "single", average = "average",
             complete = "complete", weighted = "mcquitty")
    if (method %in% names (map))
      return (unname (map [method]))
    known = c ("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
               "median", "centroid")
    if (method %in% known)
      return (method)
    stop ("HCA: unknown linkage \"", method, "\". Available: ",
          paste (sort (unique (c (names (map), known))), collapse = ", "), ".")
  }

#' Hierarchical Cluster Analysis method
#'
#' Run the HCA method for clustering.
#' @name HCA
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of cluster. If \code{NULL} (the default), it is set to the largest drop
#' in aggregation height.
#' @param method Character string defining the clustering method.
#' @param engine Which implementation builds the hierarchy: \code{\link[stats]{hclust}} (the
#' default) or \code{\link[cluster]{agnes}}. They give the same hierarchy -- same heights, to
#' machine precision, and the same cut -- but \code{hclust} is far faster on a large dataset
#' (0.3 s against 30 s on 3000 observations). Use \code{"agnes"} for the linkages it alone
#' provides, or to compare the two.
#' @param graph A logical indicating whether or not a graphic should be plotted (the
#' aggregation heights used to choose \code{k}, when \code{k} is not given).
#' @param ... Other parameters.
#' @return The cluster hierarchy (\code{hca} object).
#' @export
#' @seealso \code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}}, \code{\link{treeplot}},
#' \code{\link{predict.hca}}
#' @examples
#' require (datasets)
#' data (iris)
#' HCA (iris [, -5], k = 3, method = "ward")
HCA <-
  function (d, k = NULL, method = c ("ward", "single"), engine = c ("hclust", "agnes"),
            graph = FALSE, ...)
  {
    engine = match.arg (engine [1], c ("hclust", "agnes"))
    method = method [1]
    if (engine == "agnes")
    {
      if (!requireNamespace ("cluster", quietly = TRUE))
        stop ("HCA: engine = \"agnes\" needs the 'cluster' package. Please install it, or use ",
              "engine = \"hclust\".")
      hc = stats::as.hclust (cluster::agnes (d, method = method))
    }
    else
      hc = stats::hclust (stats::dist (d), method = hca.linkage (method))
    if (is.null (k))
    {
      heights = sort (hc$height, decreasing = TRUE)
      k = 1 + which.min (diff (heights))
      if (graph)
      {
        n = min (length (heights), 20)
        graphics::plot (1:n, heights [1:n], t = "b", col = "darkblue", xlab = "Merge",
                        ylab = "Aggregation height",
                        main = paste ("HCA: k =", k, "(largest drop in height)"))
        graphics::abline (v = k - 1, lty = 2, col = "blue")
      }
    }
    cluster = stats::cutree (hc, k)
    # The dataset travels with the result, as it does for every other clustering of the
    # package: predict.hca() needs it to place the cluster centres.
    r = c (hc, list (cluster = cluster, k = k, data = as.matrix (d)))
    class (r) = c ("hca", class (hc))
    return (r)
  }

#' @keywords internal
# Lookup table for the intern.* dispatch used by intern(). Built lazily (as a function
# rather than a static list) so it does not depend on file collation order.
intern.functions <-
  function ()
  {
    list (dunn = intern.dunn,
          interclass = intern.interclass,
          intraclass = intern.intraclass)
  }

#' Clustering evaluation through internal criteria
#'
#' Evaluation a clustering algorithm according to internal criteria.
#' @name intern
#' @param clus The extracted clusters.
#' @param d The dataset.
#' @param eval The evaluation criteria.
#' @param type Indicates whether a "global" or a "cluster"-wise evaluation should be used.
#' @return The evaluation of the clustering.
#' @export
#' @seealso \code{\link{compare}}, \code{\link{stability}}, \code{\link{intern.dunn}}, \code{\link{intern.interclass}}, \code{\link{intern.intraclass}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' intern (km$cluster, iris [, -5])
#' intern (km$cluster, iris [, -5], type = "cluster")
#' intern (km$cluster, iris [, -5], eval = c ("intraclass", "interclass"))
#' intern (km$cluster, iris [, -5], eval = c ("intraclass", "interclass"), type = "cluster")
intern <-
  function (clus, d, eval = "intraclass", type = c ("global", "cluster"))
  {
    type = match.arg (type [1], c ("global", "cluster"))
    funs = intern.functions ()
    res = sapply (eval, function (e)
    {
      fun = funs [[e]]
      if (is.null (fun))
        stop ("intern: unknown evaluation criterion '", e, "'. Available criteria: ",
              paste (names (funs), collapse = ", "))
      return (fun (clus, d, type))
    })
    if (is.vector (res))
    {
      if (type [1] == "global")
        names (res) = eval
      else
        names (res) = paste ("Cluster", sort (unique (clus)))
    }
    else
    {
      res = t (res)
      colnames (res) = paste ("Cluster", sort (unique (clus)))
      rownames (res) = eval
    }
    return (res)
  }

#' Clustering evaluation through Dunn's index
#'
#' Evaluation a clustering algorithm according to Dunn's index.
#' @name intern.dunn
#' @param clus The extracted clusters.
#' @param d The dataset.
#' @param type Indicates whether a "global" or a "cluster"-wise evaluation should be used. The
#' per-cluster values are the terms the global index is the minimum of: each cluster's distance
#' to the nearest other one, over the largest diameter of the partition.
#' @return The evaluation of the clustering.
#' @export
#' @seealso \code{\link{intern}}, \code{\link{intern.interclass}}, \code{\link{intern.intraclass}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' intern.dunn (km$cluster, iris [, -5])
#' intern.dunn (km$cluster, iris [, -5], type = "cluster")
intern.dunn <-
  function (clus, d, type = c ("global", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    type = match.arg (type [1], c ("global", "cluster"))
    clusters = sort (unique (clus))
    # A single cluster has nothing to be separated from: the inner sapply() then ran over an
    # empty set and min() returned Inf after warning about an empty argument, so the index was
    # a meaningless Inf preceded by a warning nobody could act on.
    if (length (clusters) < 2)
      stop ("intern.dunn: Dunn's index compares the distance *between* clusters with their ",
            "diameters, so it needs at least two of them; 'clus' has ", length (clusters), ".")
    # Cluster by cluster rather than on one n x n matrix: the same distances are computed,
    # but only one block of them is held at a time. On 5000 observations that is 300 Mo down
    # to a few, and the diameters need the within-cluster blocks only.
    d = as.matrix (d)
    idx = lapply (clusters, function (i) which (clus == i))
    diam = sapply (idx, function (i)
      if (length (i) < 2) 0 else max (flexclust::dist2 (d [i, , drop = FALSE], d [i, , drop = FALSE])))
    sep = matrix (Inf, nrow = length (clusters), ncol = length (clusters))
    for (a in seq_along (clusters))
      for (b in seq_len (a - 1))
      {
        closest = min (flexclust::dist2 (d [idx [[a]], , drop = FALSE], d [idx [[b]], , drop = FALSE]))
        sep [a, b] = closest
        sep [b, a] = closest
      }
    # The diagonal stays at Inf, which never wins a minimum: there are at least two clusters.
    dmin = apply (sep, 1, min)
    # Both readings divide by the *largest* diameter, so that the per-cluster values decompose
    # the global index and its minimum is the global one. Dividing each separation by its own
    # cluster's diameter made them two different statistics under one name -- and gave Inf for
    # a cluster of a single observation, whose diameter is zero.
    largest = max (diam)
    if (largest == 0)
    {
      # Every cluster reduced to a single observation: separated points with no extent are as
      # good as a clustering gets, so the limit is Inf rather than 0 / 0.
      res = rep (Inf, length (clusters))
    }
    else
      res = dmin / largest
    if (type == "global")
      return (min (res))
    names (res) = paste ("Cluster", clusters)
    return (res)
  }

#' Clustering evaluation through interclass inertia
#'
#' Evaluation a clustering algorithm according to interclass inertia.
#' @name intern.interclass
#' @param clus The extracted clusters.
#' @param d The dataset.
#' @param type Indicates whether a "global" or a "cluster"-wise evaluation should be used.
#' @return The evaluation of the clustering.
#' @export
#' @seealso \code{\link{intern}}, \code{\link{intern.dunn}}, \code{\link{intern.intraclass}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' intern.interclass (km$cluster, iris [, -5])
intern.interclass <-
  function (clus, d, type = c ("global", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    type = match.arg (type [1], c ("global", "cluster"))
    centers = apply (d, 2, function (v) tapply (v, clus, mean))
    center = matrix (apply (d, 2, mean), nrow = 1)
    res = flexclust::dist2 (center, centers)^2 * as.numeric (table (clus))
    if (type [1] == "global")
      res = sum (res)
    return (res)
  }

#' Clustering evaluation through intraclass inertia
#'
#' Evaluation a clustering algorithm according to intraclass inertia.
#' @name intern.intraclass
#' @param clus The extracted clusters.
#' @param d The dataset.
#' @param type Indicates whether a "global" or a "cluster"-wise evaluation should be used.
#' @return The evaluation of the clustering.
#' @export
#' @seealso \code{\link{intern}}, \code{\link{intern.dunn}}, \code{\link{intern.interclass}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' intern.intraclass (km$cluster, iris [, -5])
intern.intraclass <-
  function (clus, d, type = c ("global", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    type = match.arg (type [1], c ("global", "cluster"))
    centers = apply (d, 2, function (v) tapply (v, clus, mean))
    indices = sort (unique (clus))
    # tapply() returns one row per *existing* cluster, in sorted order, so the rows of
    # 'centers' are indexed by the rank of a cluster, not by its value. Using the value
    # picked the wrong row (or none at all) as soon as the ids were not exactly 1..k --
    # which is the normal case for DBSCAN, whose noise points are labelled 0.
    res = sapply (seq_along (indices), function (j)
    {
      center = matrix (centers [j, ], nrow = 1)
      return (sum (flexclust::dist2 (center, d [clus == indices [j], , drop = FALSE])^2))
    })
    if (type == "global")
      res = sum (res)
    return (res)
  }

#' @keywords internal
jaccard0 <-
  function (clus, gt)
  {
    if (!is.vector (clus))
      clus = clus$cluster
    return (sum (clus & gt) / sum (clus | gt))
  }

#' @keywords internal
jaccard1 <-
  function (clus, gt)
  {
    if (!is.vector (clus))
      clus = clus$cluster
    if (sum (clus) == 0)
      return (NA)
    return (bestoverclasses (jaccard0, clus, gt))
  }

#' @keywords internal
kappa0 <-
  function (clus, gt)
  {
    if (!is.vector (clus))
      clus = clus$cluster
    return (irr::kappa2 (cbind (clus, gt), weight = "equal")$value)
  }

#' @keywords internal
kappa1 <-
  function (clus, gt)
  {
    if (!is.vector (clus))
      clus = clus$cluster
    return (bestoverclasses (kappa0, clus, gt))
  }

#' K-means method
#'
#' Run K-means for clustering.
#'
#' The four criteria \code{criterion} offers, all computed between 2 clusters and \code{k}:
#' \describe{
#'   \item{\code{"pseudo-F"}}{the Calinski-Harabasz index, between-cluster over within-cluster
#'     variance corrected for the number of clusters. Maximised.}
#'   \item{\code{"silhouette"}}{the mean silhouette width -- how much closer each observation
#'     is to its own cluster than to the nearest other one. Maximised.}
#'   \item{\code{"gap"}}{the gap statistic: the distance between the observed within-cluster
#'     dispersion and the one expected with no cluster structure at all. The retained \code{k}
#'     is the smallest whose gap is within one standard error of the next. It is the only
#'     criterion that can answer \code{k = 1}, and much the slowest, needing \code{B}
#'     bootstrap samples.}
#'   \item{\code{"elbow"}}{the bend of the total within-cluster sum of squares. That quantity
#'     decreases with \code{k} whatever the data, so there is no optimum to take: the retained
#'     \code{k} is the point furthest from the chord joining the two ends of the curve, drawn
#'     on the graphic.}
#' }
#' The last three need the \pkg{cluster} package.
#' @name KMEANS
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of cluster.
#' @param criterion How the number of clusters is chosen: \code{"none"} (the default, use
#' \code{k} as it is), \code{"pseudo-F"}, \code{"silhouette"}, \code{"gap"} or
#' \code{"elbow"}. With any but the first, \code{k} is read as the \emph{largest} number of
#' clusters to consider. See the Details section.
#' @param B The number of bootstrap samples used by \code{criterion = "gap"}.
#' @param graph A logical indicating whether or not a graphic should be plotted (cluster number selection).
#' @param nstart Define how many random sets should be chosen.
#' @param seed A specified seed for random number generation. \emph{K}-means starts from a
#' random initialisation, so without a seed two calls on the same data give different
#' clusterings; every other clustering function of the package already had this parameter.
#' @param ... Other parameters.
#' @return The clustering (\code{kmeans} object).
#' @export
#' @seealso \code{\link[stats]{kmeans}}, \code{\link{predict.kmeans}}
#' @examples
#' require (datasets)
#' data (iris)
#' KMEANS (iris [, -5], k = 3)
#' KMEANS (iris [, -5], criterion = "pseudo-F") # With automatic detection of the nmber of clusters
KMEANS <-
  function (d, k = 9, criterion = c ("none", "pseudo-F", "silhouette", "gap", "elbow"),
            nstart = 10, B = 100, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    criterion = match.arg (criterion [1], c ("none", "pseudo-F", "silhouette", "gap", "elbow"))
    kk = k
    if (criterion != "none")
      kk = kmeans.getk (d = d, max = k, criterion = criterion, nstart = nstart, B = B,
                        graph = graph)
    res = stats::kmeans (d, kk, nstart = nstart)
    # As every other clustering of the package does.
    res$data = as.matrix (d)
    return (res)
  }

#' Estimation of the number of clusters for \emph{K}-means
#'
#' Estimate the optimal number of cluster of the \emph{K}-means clustering method.
#' @name kmeans.getk
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param max The largest number of clusters considered. Values from 2 to \code{max} are
#' evaluated (from 1, for \code{"gap"} and \code{"elbow"}, which are defined there).
#' @inheritParams KMEANS
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param nstart The number of random sets chosen for \code{\link[stats]{kmeans}} initialization.
#' @param seed A specified seed for random number generation.
#' @return The number of clusters retained by the chosen criterion.
#' @export
#' @seealso \code{\link{pseudoF}}, \code{\link{KMEANS}}, \code{\link[stats]{kmeans}},
#' \code{\link[cluster]{silhouette}}, \code{\link[cluster]{clusGap}}
#' @references Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number of
#' clusters in a data set via the gap statistic. \emph{Journal of the Royal Statistical
#' Society: Series B}, 63(2), 411-423.
#' @examples
#' require (datasets)
#' data (iris)
#' kmeans.getk (iris [, -5])
#' kmeans.getk (iris [, -5], criterion = "silhouette")
#' kmeans.getk (iris [, -5], criterion = "elbow")
#' \donttest{
#' # The gap statistic resamples, so it is much slower than the other three.
#' kmeans.getk (iris [, -5], criterion = "gap", B = 20, seed = 0)
#' }
kmeans.getk <-
  function (d, max = 9, criterion = c ("pseudo-F", "silhouette", "gap", "elbow"),
            nstart = 10, B = 100, graph = FALSE, seed = NULL)
  {
    setseed (seed)
    # 'max' below 2 would make the loops below run backwards (2:1 is c (2, 1)): they would fit
    # a two-cluster partition, then a one-cluster one, and return k = 2 regardless of 'max'.
    if (max < 2)
      stop ("kmeans.getk: 'max' must be at least 2 to compare partitions, but is ", max, ".")
    criterion = match.arg (criterion [1], c ("pseudo-F", "silhouette", "gap", "elbow"))
    fit = function (k) stats::kmeans (d, k, nstart = nstart)
    k = NA
    sizes = 2:max
    measure = NULL
    measure2 = NULL
    criterion2 = NULL
    extra = NULL
    if (criterion == "pseudo-F")
    {
      # Calinski-Harabasz: between-cluster variance over within-cluster variance, corrected
      # for the number of clusters. Maximised.
      measure2 = vector ("numeric", max - 1)
      criterion2 = as.expression (substitute (R^2))
      measure = vector ("numeric", max - 1)
      for (i in sizes)
      {
        km = fit (i)
        measure [i - 1] = pseudoF (km)
        measure2 [i - 1] = km$betweenss / km$totss
      }
      k = sizes [which.max (measure)]
    }
    else if (criterion == "silhouette")
    {
      # Mean silhouette width: for each observation, how much closer it is to its own cluster
      # than to the nearest other one, between -1 and 1. Maximised.
      if (!requireNamespace ("cluster", quietly = TRUE))
        stop ("kmeans.getk: criterion = \"silhouette\" needs the 'cluster' package. ",
              "Please install it, or use another criterion.")
      dis = stats::dist (d)
      measure = sapply (sizes, function (i)
        mean (cluster::silhouette (fit (i)$cluster, dis) [, "sil_width"]))
      k = sizes [which.max (measure)]
    }
    else if (criterion == "gap")
    {
      # Tibshirani, Walther & Hastie (2001): the gap between the observed within-cluster
      # dispersion and the one expected under a null model with no cluster structure, both on
      # a log scale. The retained k is the smallest one whose gap is within one standard error
      # of the following one -- not simply the maximum.
      if (!requireNamespace ("cluster", quietly = TRUE))
        stop ("kmeans.getk: criterion = \"gap\" needs the 'cluster' package. ",
              "Please install it, or use another criterion.")
      g = cluster::clusGap (as.matrix (d), FUNcluster = stats::kmeans, K.max = max, B = B,
                            nstart = nstart, verbose = FALSE)
      sizes = 1:max
      measure = g$Tab [, "gap"]
      extra = g$Tab [, "SE.sim"]
      k = cluster::maxSE (measure, extra, method = "firstSEmax")
    }
    else
    {
      # Elbow: the total within-cluster sum of squares always decreases with k, so there is no
      # optimum to take -- one looks for the bend. knee.index() picks the point furthest from
      # the chord joining the two ends of the curve, i.e. the bend one would point at by eye.
      sizes = 1:max
      measure = sapply (sizes, function (i) fit (i)$tot.withinss)
      k = sizes [knee.index (measure)]
    }
    if (graph & !is.na (k))
    {
      ylab = c ("pseudo-F" = "pseudo-F", silhouette = "Mean silhouette width",
                gap = "Gap statistic", elbow = "Total within-cluster sum of squares") [criterion]
      if (!is.null (measure2))
      {
        opar1 = graphics::par (mar = c(5, 4, 4, 5) + .1)
        on.exit (graphics::par (opar1), add = TRUE)
      }
      ylim = range (measure)
      if (!is.null (extra))
        ylim = range (c (measure - extra, measure + extra))
      graphics::plot (sizes, measure, t = "b", xlab = "Number of clusters", ylab = ylab,
                      ylim = ylim)
      if (!is.null (extra))
        graphics::arrows (sizes, measure - extra, sizes, measure + extra, length = .03,
                          angle = 90, code = 3, col = "grey40")
      if (criterion == "elbow")
        # The chord the bend is measured against, so that the construction is visible.
        graphics::segments (sizes [1], measure [1], sizes [length (sizes)],
                            measure [length (measure)], lty = 2, col = "grey40")
      if (!is.null (measure2))
      {
        opar2 = graphics::par (new = TRUE)
        on.exit (graphics::par (opar2), add = TRUE)
        graphics::plot (2:max, measure2, xaxt="n", yaxt="n", xlab="", ylab="", t= "b", lty = 2, pch = 2)
        graphics::axis (4)
        graphics::mtext (criterion2, side = 4, line = 3)
        graphics::legend ("bottomright", c (criterion, criterion2), lty = 1:2)
      }
      graphics::abline (v = k, lty = 3)
    }
    return (k)
  }

#' MeanShift method
#'
#' Run MeanShift for clustering.
#' @name MEANSHIFT
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param mskernel A string indicating the kernel associated with the kernel density estimate that the mean shift is optimizing over.
#' @param bandwidth Used in the kernel density estimate for steepest ascent classification.
#' @param alpha A scalar tuning parameter for normal kernels.
#' @param iterations The number of iterations to perform mean shift.
#' @param epsilon A scalar used to determine when to terminate the iteration of an individual query point.
#' @param epsilonCluster A scalar used to determine the minimum distance between distinct clusters.
#' @param seed A specified seed for random number generation. The MeanShift algorithm itself is deterministic
#' given its parameters, but the seed is provided for consistency with the rest of the package's API.
#' @param ... Other parameters.
#' @return The clustering (\code{meanshift} object).
#' @export
#' @seealso \code{\link[meanShiftR]{meanShift}}, \code{\link{predict.meanshift}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' MEANSHIFT (iris [, -5], bandwidth = .75)
#' }
MEANSHIFT <-
  function (d, mskernel = "NORMAL", bandwidth = rep (1, ncol (d)), alpha = 0, iterations = 10, epsilon = 1e-08, epsilonCluster = 1e-04, seed = NULL, ...)
  {
    setseed (seed)
    dd = as.matrix (d)
    if (length (bandwidth) == 1)
      bandwidth = rep (bandwidth, ncol (d))
    res = meanShiftR::meanShift (dd, kernelType = mskernel, bandwidth = bandwidth, alpha = alpha, iterations = iterations, epsilon = epsilon, epsilonCluster = epsilonCluster)
    names (res) [1] = "cluster"
    res [[1]] = as.vector (res [[1]])
    res = c (res,
             data = list (dd),
             kernel = list (mskernel),
             bandwidth = list (bandwidth),
             alpha = list (alpha),
             iterations = list (iterations),
             epsilon = list (epsilon),
             epsilonCluster = list (epsilonCluster))
    class (res) = "meanshift"
    return (res)
  }

#' Plot function for som-class
#'
#' Plot Kohonen's self-organizing maps.
#' @name plot.som
#' @param x The Kohonen's map (object of class \code{\link{som-class}}).
#' @param type The type of plot.
#' @param col Color of the data points
#' @param labels A \code{vector} of character strings to be printed instead of points in the plot.
#' @param ... Other parameters.
#' @export
#' @method plot som
#' @seealso \code{\link{SOM}}, \code{\link{som-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' som = SOM (iris [, -5], xdim = 5, ydim = 5, post = "ward", k = 3)
#' plot (som) # Scatter plot (default)
#' plot (som, type = "mapping") # Kohonen map
plot.som <-
  function (x, type = c ("scatter", "mapping"), col = NULL, labels = FALSE, ...)
  {
    if (type [1] == "scatter")
    {
      d = x$data
      # kohonen stores the codebook vectors in a *list*, one element per data layer, so the
      # matrix is codes [[1]] and not codes.
      centers.coord = x$som$codes [[1]]
      xlab = colnames (d) [1]
      ylab = colnames (d) [2]
      col = x$cluster + 1
      if (ncol (x$data) > 2)
      {
        proj = pca.project2d (d)
        centers.coord = proj$project (x$som$codes [[1]])
        d = proj$coord
        xlab = proj$xlab
        ylab = proj$ylab
      }
      if ((!labels) | (is.null (row.names (d))))
        graphics::plot (d, col = col, xlab = xlab, ylab = ylab, asp = 1)
      else
      {
        graphics::plot (d, col = 0, xlab = xlab, ylab = ylab, asp = 1)
        graphics::text (d, row.names (d), col = col)
      }
      pts = x$som$grid$pts
      adj = which (as.matrix (stats::dist (pts [, 1])) ==  1 & as.matrix (stats::dist (pts [, 2])) ==  0, arr.ind = TRUE)
      # which (arr.ind = TRUE) returns both (i, j) and (j, i): keep one of each pair, or every
      # segment gets drawn twice.
      adj = adj [adj [, 1] < adj [, 2], , drop = FALSE]
      graphics::segments (centers.coord [adj [, 1], 1], centers.coord [adj [, 1], 2],
                          centers.coord [adj [, 2], 1], centers.coord [adj [, 2], 2],
                          lty = 2, col = "darkgrey")
      adj = which (as.matrix (stats::dist (pts [, 1])) ==  0 & as.matrix (stats::dist (pts [, 2])) ==  1, arr.ind = TRUE)
      adj = adj [adj [, 1] < adj [, 2], , drop = FALSE]
      graphics::segments (centers.coord [adj [, 1], 1], centers.coord [adj [, 1], 2],
                          centers.coord [adj [, 2], 1], centers.coord [adj [, 2], 2],
                          lty = 2, col = "darkgrey")
      graphics::points (centers.coord, col = x$nodes + 1, pch = 19)
    }
    else if (type [1] == "mapping")
    {
      if (!is.null (col))
      {
        d = x$data
        pcol = col
        tmp = table (col, factor (x$som$unit.classif, levels = min (x$som$unit.classif):max(x$som$unit.classif)))
        bgcol = apply (tmp, 2, which.max) + 1
        empty = apply (tmp, 2, max) == 0
        bgcol [empty] = col [apply (flexclust::dist2 (d, x$som$codes [[1]] [empty, , drop = FALSE]), 2, which.min)]
        bgcol = grDevices::rgb (t ((grDevices::col2rgb (bgcol) * 2 + 255) / 3), maxColorValue = 255)
      }
      else
      {
        pcol = 1
        bgcol = 0
        if (length (unique (x$nodes)) != length (x$nodes))
        {
          bgcol = grDevices::rgb (t ((grDevices::col2rgb (x$nodes + 1) * 2 + 255) / 3), maxColorValue = 255)
          pcol = x$cluster + 1
        }
      }
      lab = NULL
      if (labels)
        lab = row.names (x$data)
      # x is the fdm2id 'som' object (som / nodes / cluster / data); the unit assignments live
      # on the wrapped kohonen map, so x$unit.classif was always NULL here.
      graphics::plot (x$som, type = "mapping", classif = x$som$unit.classif,
                      pch = 19, col = pcol, bgcol = bgcol, main = "", labels = lab)
    }
    else
      message (paste (type, ": Unknown"))
  }

#' @keywords internal
plotclus.extract.kmeans <-
  function (clustering, centers, k)
  {
    list (clusters = clustering$cluster, centres = if (centers) clustering$centers else NULL, k = k)
  }

#' @keywords internal
plotclus.extract.em <-
  function (clustering, centers, k)
  {
    list (clusters = apply (clustering$z, 1, which.max),
         centres = if (centers) t (clustering$parameters$mean) else NULL, k = k)
  }

#' @keywords internal
# Shared extractor for the clustering methods that only expose a plain $cluster vector and
# no notion of centers (DBSCAN, SOM, MeanShift, Spectral).
plotclus.extract.simple <-
  function (clustering, centers, k)
  {
    list (clusters = clustering$cluster, centres = NULL, k = k)
  }

#' @keywords internal
plotclus.extract.hca <-
  function (clustering, centers, k)
  {
    if (is.null (k))
    {
      if (is.null (clustering$cluster))
        k = 1 + which.min (diff (sort (clustering$height, decreasing = TRUE)))
      else
        k = length (unique (clustering$cluster))
    }
    list (clusters = stats::cutree (clustering, k), centres = NULL, k = k)
  }

#' @keywords internal
# Lookup table mapping a clustering-method S3/S4 class to the function that extracts its
# cluster assignments (and, when available and requested, cluster centers) for plotclus().
# Built lazily (as a function rather than a static list) so it does not depend on file
# collation order. Order matters: it mirrors the priority of the original if/else-if chain
# (only the first matching class in this list is used).
plotclus.extractors <-
  function ()
  {
    list (kmeans = plotclus.extract.kmeans,
          pam = plotclus.extract.kmeans,
          em = plotclus.extract.em,
          dbs = plotclus.extract.simple,
          som = plotclus.extract.simple,
          meanshift = plotclus.extract.simple,
          spectral = plotclus.extract.simple,
          hca = plotclus.extract.hca)
  }

#' Generic Plot Method for Clustering
#'
#' Plot a clustering according to various parameters
#' @name plotclus
#' @param clustering The clustering to be plotted.
#' @param d The dataset (\code{matrix} or \code{data.frame}), mandatory for some of the graphics.
#' @param type The type of plot.
#' @param centers Indicates whether or not cluster centers should be plotted (used only in scatter plots).
#' @param k Number of clusters (used only for hierarchical methods). If not specified an "optimal" value is determined.
#' @param tailsize Number of clusters showned (used only for height plots).
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{treeplot}}, \code{\link{scatterplot}}, \code{\link{plot.som}}, \code{\link{boxclus}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' ward = HCA (iris [, -5], k = 3, method = "ward")
#' plotclus (ward, iris [, -5], type = "scatter") # Scatter plot
#' plotclus (ward, iris [, -5], type = "boxplot") # Boxplot
#' plotclus (ward, iris [, -5], type = "tree") # Dendrogram
#' plotclus (ward, iris [, -5], type = "height") # Distances between merging clusters
#' som = SOM (iris [, -5], xdim = 5, ydim = 5, post = "ward", k = 3)
#' plotclus (som, iris [, -5], type = "scatter") # Scatter plot for SOM
#' plotclus (som, iris [, -5], type = "mapping") # Kohonen map
#' }
plotclus <-
  function (clustering,
            d = NULL,
            type = c ("scatter", "boxplot", "tree", "height", "mapping", "words"),
            centers = FALSE,
            k = NULL,
            tailsize = 9,
            ...)
  {
    type = match.arg (type [1], c ("scatter", "boxplot", "tree", "height", "mapping", "words"))
    method = class (clustering)
    extractors = plotclus.extractors ()
    matched = intersect (names (extractors), method)
    clusters = NULL
    centres = NULL
    if (length (matched) > 0)
    {
      extracted = extractors [[matched [1]]] (clustering, centers, k)
      clusters = extracted$clusters
      centres = extracted$centres
      k = extracted$k
    }
    ishca = "hca" %in% method
    issom = "som" %in% method
    if ((type [1] == "tree") & ishca)
      treeplot (clustering, k = k, ...)
    else if (type [1] == "height" & ishca)
      graphics::barplot (utils::tail (sort (clustering$height), n = tailsize), names.arg = tailsize:1,
                         xlab = "Number of clusters", ylab = "Height", main = "", sub = "")
    else if ((type [1] == "scatter") & (!is.vector (d)) & issom)
      plot.som (clustering, type = type, ...)
    else if ((type [1] == "mapping") & issom)
      plot.som (clustering, type = type, ...)
    else if (type [1] == "scatter")
      scatterplot (d, clusters, centres, ...)
    else if (type [1] == "boxplot")
      boxclus (d, clusters, ...)
    else if (type [1] == "words")
      plotcloud (d, clusters, ...)
    else
      stop ("plotclus: type = \"", type, "\" is not available for a clustering of class \"",
            method [1], "\". Dendrograms and aggregation heights (\"tree\", \"height\") need a ",
            "hierarchical clustering, and \"mapping\" a self-organising map.")
  }

#' Predict function for DBSCAN
#'
#' Return the closest DBSCAN cluster for a new dataset.
#' @name predict.dbs
#' @param object The classification model (of class \code{\link{dbs-class}}, created by \code{\link{DBSCAN}}).
#' @param newdata A new dataset (a \code{data.frame}), with same variables as the learning dataset.
#' @param ... Other parameters.
#' @export
#' @method predict dbs
#' @seealso \code{\link{DBSCAN}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = DBSCAN (d$train.x, minpts = 5, eps = 0.65)
#' predict (model, d$test.x)
predict.dbs <-
  function (object, newdata, ...)
  {
    select = object$isseed | (object$cluster == 0)
    object$cluster [select][apply (flexclust::dist2 (object$data [select, ], newdata), 2, which.min)]
  }

#' Predict function for EM
#'
#' Return the closest EM cluster for a new dataset.
#' @name predict.em
#' @param object The classification model (of class \code{\link{em-class}}, created by \code{\link{EM}}).
#' @param newdata A new dataset (a \code{data.frame}), with same variables as the learning dataset.
#' @param ... Other parameters.
#' @export
#' @method predict em
#' @seealso \code{\link{EM}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = EM (d$train.x, 3)
#' predict (model, d$test.x)
predict.em <-
  function (object, newdata, ...)
  {
    apply (mclust::estep (data = newdata, modelName = object$modelName, parameters = object$parameters)$z, 1, which.max)
  }

#' Predict function for K-means
#'
#' Return the closest K-means cluster for a new dataset.
#' @name predict.kmeans
#' @param object The classification model (created by \code{\link{KMEANS}}).
#' @param newdata A new dataset (a \code{data.frame}), with same variables as the learning dataset.
#' @param ... Other parameters.
#' @export
#' @method predict kmeans
#' @seealso \code{\link{KMEANS}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = KMEANS (d$train.x, k = 3)
#' predict (model, d$test.x)
predict.kmeans <-
  function (object, newdata, ...)
  {
    apply (flexclust::dist2 (object$centers, newdata), 2, which.min)
  }


#' Predict function for MeanShift
#'
#' Return the closest MeanShift cluster for a new dataset.
#' @name predict.meanshift
#' @param object The classification model (created by \code{\link{MEANSHIFT}}).
#' @param newdata A new dataset (a \code{data.frame}), with same variables as the learning dataset.
#' @param ... Other parameters.
#' @export
#' @method predict meanshift
#' @seealso \code{\link{MEANSHIFT}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = MEANSHIFT (d$train.x, bandwidth = .75)
#' predict (model, d$test.x)
#' }
predict.meanshift <-
  function (object, newdata, ...)
  {
    res = meanShiftR::meanShift (queryData = as.matrix (newdata),
                                 trainData = object$data,
                                 kernelType = object$kernel,
                                 bandwidth = object$bandwidth,
                                 alpha = object$alpha,
                                 iterations = object$iterations,
                                 epsilon = object$epsilon,
                                 epsilonCluster = object$epsilonCluster)
    mmodel = apply (object$value, 2, function (v) tapply (v, object$cluster, mean))
    mpred = apply (res$value, 2, function (v) tapply (v, res$assignment, mean))
    conv = as.vector (apply (flexclust::dist2 (mmodel, mpred), 2, which.min))
    return (conv [as.vector (res$assignment)])
  }

#' Pseudo-F
#'
#' Compute the pseudo-F of a clustering result obtained by the \emph{K}-means method.
#' @name pseudoF
#' @param clustering The clustering result (obtained by the function \code{\link[stats]{kmeans}}).
#' @return The pseudo-F of the clustering result.
#' @export
#' @seealso \code{\link{kmeans.getk}}, \code{\link{KMEANS}}, \code{\link[stats]{kmeans}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' pseudoF (km)
pseudoF <-
  function (clustering)
  {
    r2 = clustering$betweenss / clustering$totss
    k = length (clustering$size)
    n = length (clustering$cluster)
    f = (r2 / (k-1)) / ((1-r2) / (n-k))
    return (f)
  }

#' Clustering Scatter Plots
#'
#' Produce a scatter plot for clustering results. If the dataset has more than two dimensions, the scatter plot will show the two first PCA axes.
#' @name scatterplot
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param clusters Cluster labels of the training set: a numeric \code{vector} (0 marking the
#' observations a density method left as noise), or a \code{factor}/character vector, whose
#' levels are then used to label the legend.
#' @param centers Coordinates of the cluster centers.
#' @param labels Indicates whether or not labels (row names) should be showned on the plot.
#' @param ellipses Indicates whether or not ellipses should be drawned around clusters.
#' @param legend Indicates where the legend is placed on the graphics.
#' @param ... Other parameters.
#' @export
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' scatterplot (iris [, -5], km$cluster)
scatterplot <-
  function (d, clusters, centers = NULL, labels = FALSE, ellipses = FALSE, legend = c ("auto1", "auto2"), ...)
  {
    # 'clusters' is documented as a vector *or a factor*, but the whole body treats it as
    # numeric -- min(), max() and 1 + clusters. Categorical input (the ground truth of a
    # dataset, the output of a classifier: the natural things to colour a scatter plot by) is
    # converted to the integer codes the rest of the function expects, its labels being kept
    # for the legend. Shared with boxclus(), which draws the same clusters.
    coded = cluster.codes (clusters)
    clusters = coded$codes
    clusternames = coded$names
    dd = NULL
    kmin = 1 + min (clusters)
    kmax = 1 + max (clusters)
    col = 1 + clusters
    asp = TRUE
    if (is.vector (d))
    {
      d = cbind (Index = 1:(length (d)), Data = d)
      asp = FALSE
    }
    if (ncol (d) == 2)
    {
      dd = d
      xlab = colnames (d) [1]
      ylab = colnames (d) [2]
    }
    else
    {
      proj = pca.project2d (d)
      dd = proj$coord
      if (!is.null (centers))
        centers = proj$project (centers)
      xlab = proj$xlab
      ylab = proj$ylab
    }
    if ((!labels) | (is.null (row.names (dd))))
    {
      if (asp)
        graphics::plot (dd, asp = 1, col = col, xlab = xlab, ylab = ylab)
      else
        graphics::plot (dd, col = col, xlab = xlab, ylab = ylab)
    }
    else
    {
      if (asp)
        graphics::plot (dd, asp = 1, col = 0, xlab = xlab, ylab = ylab)
      else
        graphics::plot (dd, col = 0, xlab = xlab, ylab = ylab)
      graphics::text (dd, row.names (dd), col = col)
    }
    if (!is.null (centers))
      graphics::points (centers [, 1], centers [, 2], pch = 19, col = kmin:kmax)
    if (ellipses)
    {
      # mclust::unmap() numbers the mixture components 1..k in the order of the sorted
      # cluster values, so sigma [,, j] and the rows of 'centers' are indexed by the *rank*
      # of a cluster, not by its value. Looping over the values themselves went out of
      # bounds as soon as the ids had a gap (clusters 1 and 3, say). 'ecol' rather than
      # 'col', which still holds the per-observation colours the legend below needs.
      keep = clusters != 0
      lev = sort (unique (clusters [keep]))
      z = mclust::unmap (clusters [keep])
      p = mclust::mstep (data = dd [keep, ], modelName = "VVV", z = z)
      if (is.null (centers))
        centers = t (p$parameters$mean)
      for (j in seq_along (lev))
      {
        ecol = lev [j] + 1
        # car::ellipse (center, shape, radius, ...) -- the centre first, then the covariance
        # matrix. draw = FALSE returns the coordinates instead of plotting them itself.
        e = car::ellipse (center = centers [j, 1:2],
                          shape = p$parameters$variance$sigma [,, j],
                          radius = sqrt (stats::qchisq (0.95, 2)),
                          draw = FALSE, segments = 1000)
        graphics::polygon (e, border = ecol, lty = 2)
        eig = eigen (p$parameters$variance$sigma [,, j])
        seg = sweep (eig$vectors, 1, sqrt (eig$values), FUN = "*")
        seg = seg * sqrt (stats::qchisq (0.95, 2))
        graphics::segments (seg [1, 1] + centers [j, 1], -seg [1, 2] + centers [j, 2],
                            -seg [1, 1] + centers [j, 1], seg [1, 2] + centers [j, 2],
                            col = ecol, lty = 3)
        graphics::segments (seg [2, 1] + centers [j, 1], -seg [2, 2] + centers [j, 2],
                            -seg [2, 1] + centers [j, 1], seg [2, 2] + centers [j, 2],
                            col = ecol, lty = 3)
      }
    }
    if (legend [1] == "auto1")
    {
      coord = graphics::par ("usr")
      pos = rbind (coord [c (1, 4)], coord [c (2, 4)], coord [c (1, 3)], coord [c (2, 3)])
      legend = c ("topleft", "topright", "bottomleft", "bottomright") [which.max (apply (flexclust::dist2 (pos, d [, 1:2]), 1, min))]
    }
    else if (legend [1] == "auto2")
    {
      coord = graphics::par ("usr")
      cx = mean (coord [1:2])
      cy = mean (coord [3:4])
      left = d [, 1] < cx
      right = d [, 1] > cx
      bottom = d [, 2] < cy
      top = d [, 2] > cy

      count = c (sum (top & left), sum (top & right), sum (bottom & left), sum (bottom & right))
      legend = c ("topleft", "topright", "bottomleft", "bottomright") [which.min (count)]
    }
    entries = cluster.legend (clusters, clusternames)
    graphics::legend (x = legend, legend = entries$labels, col = entries$col, pch = 1)
  }

#' Self-Organizing Maps clustering method
#'
#' Run the SOM algorithm for clustering.
#' @name SOM
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param xdim,ydim The dimensions of the grid.
#' @param rlen The number of iterations.
#' @param post The post-treatement method: \code{"none"} (None), \code{"single"} (Single link) or \code{"ward"} (Ward clustering).
#' @param k The number of cluster (only used if \code{post} is different from \code{"none"}).
#' @param seed A specified seed for random number generation (codebook initialization).
#' @return The fitted Kohonen's map as an object of class \code{som}.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{plot.som}}, \code{\link{som-class}}, \code{\link[kohonen]{som}}
#' @examples
#' require (datasets)
#' data (iris)
#' SOM (iris [, -5], xdim = 5, ydim = 5, post = "ward", k = 3)
SOM <-
  function (d, xdim = floor (sqrt (nrow (d))), ydim = floor (sqrt (nrow (d))), rlen = 10000, post = c ("none", "single", "ward"), k = NULL, seed = NULL, ...)
  {
    setseed (seed)
    grid = kohonen::somgrid (xdim, ydim, topo  = "rectangular", toroidal = FALSE)
    map = kohonen::som (as.matrix (d), grid = grid, rlen = rlen)
    nodes = rep (1, xdim * ydim)
    cluster = rep (1, nrow (d))
    if (post [1] != "none")
    {
      hc = cluster::agnes (map$codes [[1]], method = post [1])
      if (is.null (k))
        k = 1 + which.min (diff (sort (hc$height, decreasing = TRUE)))
      nodes = stats::cutree (hc, k)
      cluster = nodes [map$unit.classif]
    }
    else
    {
      cluster = map$unit.classif
      nodes = 1:(length (nodes))
    }
    r = list (som = map, nodes = nodes, cluster = cluster, data = as.matrix (d))
    class (r) = "som"
    return (r)
  }

#' @keywords internal
# The k leading eigenvectors of a symmetric matrix.
#
# eigen() computes all n of them, which is the whole cost of a spectral clustering: 10.2 s of
# the 10.5 s it takes on 2000 observations, and it grows as n^3. RSpectra::eigs_sym() computes
# the k that are wanted and nothing else -- 0.45 s on the same matrix -- and returns the same
# subspace, hence the same clustering (checked against eigen() by the test suite). It is only
# a Suggests, so eigen() remains the fallback.
#
# symmetric = TRUE is not optional in that fallback: 'l' is symmetric by construction, but
# rounding can hide it from eigen(), which then uses the general algorithm -- complex
# eigenvalues, ordered by modulus instead of by value, and a projection built on wrong axes.
spectral.eigenvectors <-
  function (l, k)
  {
    if (requireNamespace ("RSpectra", quietly = TRUE) && (k < nrow (l) - 1))
    {
      # The iterative solver can stop before it has converged on 'k' eigenvectors -- on an
      # affinity matrix close to block-diagonal, which is exactly what a bootstrap resample
      # of a well-separated dataset produces. It then returns fewer columns than asked for,
      # and everything downstream is indexed on 'k'. Fall back to the dense decomposition.
      x = suppressWarnings (RSpectra::eigs_sym (l, k)$vectors)
      if (!is.null (x) && (ncol (x) == k))
        return (x)
    }
    return (eigen (l, symmetric = TRUE)$vectors [, 1:k, drop = FALSE])
  }

#' Spectral clustering method
#'
#' Run a Spectral clustering algorithm.
#' @name SPECTRAL
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of cluster.
#' @param sigma Width of the gaussian used to build the affinity matrix.
#' @param graph A logical indicating whether or not a graphic should be plotted (projection on the spectral space of the affinity matrix).
#' @param seed A specified seed for random number generation (final k-means step).
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{spectral-class}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' SPECTRAL (iris [, -5], k = 3)
#' }
SPECTRAL <-
  function (d, k, sigma = 1, graph = FALSE, seed = NULL, ...)
  {
    setseed (seed)
    a = exp (-flexclust::dist2 (d, d)^2 / (2 * sigma * sigma))
    diag (a) = 0
    degree = rowSums (a)
    # With a small 'sigma' every affinity of a distant observation underflows to zero, its
    # degree with them, and the normalisation below divides by it.
    if (any (degree <= 0))
      stop ("SPECTRAL: ", sum (degree <= 0), " observation(s) have no affinity at all with the ",
            "rest of the sample, so the graph cannot be normalised. 'sigma' (", sigma,
            ") is too small for the scale of these data -- try a larger value.")
    p = 1 / sqrt (degree)
    # The normalised affinity matrix, as diag (p) %*% a %*% diag (p) -- but written as a
    # scaling of the rows and the columns. The two matrix products it replaces were O(n^3) on
    # a matrix known to be diagonal: 22 s of the 22.4 s a spectral clustering of 2000
    # observations took.
    l = a * outer (p, p)
    x = spectral.eigenvectors (l, k)
    proj = sweep (x, 1, sqrt (rowSums (x * x)), "/")
    colnames (proj) = paste ("Comp.", 1:k)
    rownames (proj) = rownames (d)
    km = stats::kmeans (proj, centers = k, nstart = 100)
    cluster = km$cluster
    if (graph)
      plotdata (d = proj, k = factor (paste ("Cluster", cluster)))
    res = list (cluster = cluster, proj = proj, centers = km$centers, data = as.matrix (d))
    class (res) = "spectral"
    return (res)
  }

#' Predict function for Spectral clustering
#'
#' Return the closest Spectral clustering cluster for a new dataset. New instances are assigned to
#' the cluster of their nearest neighbour in the original (training) space, since the spectral
#' projection cannot be directly extended to unseen data without recomputing the affinity matrix.
#' @name predict.spectral
#' @param object The clustering model (of class \code{\link{spectral-class}}, created by \code{\link{SPECTRAL}}).
#' @param newdata A new dataset (a \code{data.frame}), with same variables as the learning dataset.
#' @param ... Other parameters.
#' @export
#' @method predict spectral
#' @seealso \code{\link{SPECTRAL}}, \code{\link{spectral-class}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = SPECTRAL (d$train.x, k = 3)
#' predict (model, d$test.x)
#' }
predict.spectral <-
  function (object, newdata, ...)
  {
    object$cluster [apply (flexclust::dist2 (object$data, newdata), 2, which.min)]
  }

#' Clustering evaluation through stability
#'
#' Evaluation a clustering algorithm according to stability, through a bootstrap procedure.
#' @name stability
#' @param clusteringmethods The clustering methods to be evaluated.
#' @param d The dataset.
#' @param originals The original clustering.
#' @param eval The evaluation criteria.
#' @param type The comparison method.
#' @param nsampling The number of bootstrap runs.
#' @param seed A specified seed for random number generation (useful for testing different method with the same bootstap samplings).
#' @param names Method names.
#' @param graph Indicates wether or not a graphic is potted for each sample.
#' @param ... Parameters to be passed to the clustering algorithms.
#' @return The evaluation of the clustering algorithm(s) (numeric values).
#' @export
#' @seealso \code{\link{compare}}, \code{\link{intern}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' stability (KMEANS, iris [, -5], seed = 0, k = 3)
#' stability (KMEANS, iris [, -5], seed = 0, k = 3, eval = c ("jaccard", "accuracy"), type = "global")
#' stability (KMEANS, iris [, -5], seed = 0, k = 3, type = "cluster")
#' stability (KMEANS, iris [, -5], seed = 0, k = 3, eval = c ("jaccard", "accuracy"), type = "cluster")
#' stability (c (KMEANS, HCA), iris [, -5], seed = 0, k = 3)
#' stability (c (KMEANS, HCA), iris [, -5], seed = 0, k = 3,
#' eval = c ("jaccard", "accuracy"), type = "global")
#' stability (c (KMEANS, HCA), iris [, -5], seed = 0, k = 3, type = "cluster")
#' stability (c (KMEANS, HCA), iris [, -5], seed = 0, k = 3,
#' eval = c ("jaccard", "accuracy"), type = "cluster")
#' stability (KMEANS, iris [, -5], originals = KMEANS (iris [, -5], k = 3)$cluster, seed = 0, k = 3)
#' stability (KMEANS, iris [, -5], originals = KMEANS (iris [, -5], k = 3), seed = 0, k = 3)
#' }
stability <-
  function (clusteringmethods, d, originals = NULL, eval = "jaccard", type = c ("cluster", "global"), nsampling = 10, seed = NULL, names = NULL, graph = FALSE, ...)
  {
    comp = ifelse (type [1] == "cluster", "cluster", "max")
    methodNames = names
    if (is.character (clusteringmethods))
    {
      methodNames = clusteringmethods
      clusteringmethods = sapply (clusteringmethods, get)
    }
    else
    {
      if (is.null (names))
      {
        # See performance(): a list of methods held in a variable carries no names, and the
        # result used to come back unnamed rather than numbered.
        methodNames = as.character (match.call ()$clusteringmethods)
        if (length (methodNames) > 1)
          methodNames = methodNames [-1]
        if (length (methodNames) != length (clusteringmethods))
          methodNames = paste ("Method", seq_along (clusteringmethods))
      }
    }
    clusteringmethods = c (clusteringmethods)
    if (is.null (originals))
    {
      setseed (seed + 1)
      originals = sapply (clusteringmethods, function (clus) {clus (d, graph = graph, ...)$cluster})
      originals = split (originals, rep (1:ncol (originals), each = nrow (originals)))
    }
    else
      originals = list (originals)
    nb = length (clusteringmethods)
    if (length (originals) != length (clusteringmethods))
      message ("Unsuitable number of methods")
    else
    {
      indices = 1:nb
      res = lapply (indices, function (i)
      {
        clus = clusteringmethods [[i]]
        original = originals [[i]]
        if (!is.vector (original))
          original = original$cluster
        clusternames = paste ("Cluster", sort (unique (original)))
        setseed (seed)
        s = matrix (sample (nrow (d), nsampling * nrow (d), replace = TRUE), ncol = nsampling)
        # A list, not apply()'s guesswork: a bootstrap sample need not contain every cluster of
        # the reference partition, and a method that finds its own number of clusters -- DBSCAN,
        # mean shift -- does not find the same number twice. apply() then returned a *list*
        # instead of a matrix, is.vector() is TRUE for a list, and the run ended on
        # mean (list) = NA with "argument is not numeric or logical".
        runs = lapply (seq_len (nsampling), function (j)
        {
          v = s [, j]
          clusters = clus (d [v, ], graph = graph, ...)
          if (graph)
            plotclus (clusters, d [v, ])
          return (compare (original [unique (v)], clusters$cluster [!duplicated (v)],
                           eval = eval, comp = comp))
        })
        if (comp != "cluster")
        {
          # One value per criterion and per run.
          res = rowMeans (matrix (unlist (runs), nrow = length (eval)), na.rm = TRUE)
          names (res) = eval
          return (res)
        }
        # One value per cluster *of the reference partition* and per criterion. Each run is put
        # back on those clusters -- the ones it did not see stay NA -- so that the average is
        # taken cluster by cluster and not by position.
        res = sapply (eval, function (e)
        {
          byrun = sapply (runs, function (r)
          {
            values = if (is.matrix (r)) r [e, ] else r
            return (values [match (clusternames, names (values))])
          })
          return (rowMeans (matrix (byrun, nrow = length (clusternames)), na.rm = TRUE))
        })
        res = matrix (res, nrow = length (clusternames),
                      dimnames = list (clusternames, eval))
        return (res)
      })
      if (length (res) == 1)
        res = res [[1]]
      else
        names (res) = methodNames
      return (res)
    }
  }

#' Dendrogram Plots
#'
#' Draws a dendrogram.
#' @name treeplot
#' @param clustering The dendrogram to be plotted (result of \code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}} or \code{\link{HCA}}).
#' @param labels Indicates whether or not labels (row names) should be showned on the plot.
#' @param k Number of clusters. If not specified an "optimal" value is determined.
#' @param split Indicates wheather or not the clusters should be highlighted in the graphics.
#' @param horiz Indicates if the dendrogram should be drawn horizontally or not.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link[stats]{dendrogram}}, \code{\link{HCA}}, \code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}}
#' @examples
#' require (datasets)
#' data (iris)
#' hca = HCA (iris [, -5], k = 3, method = "ward")
#' treeplot (hca)
treeplot <-
  function (clustering,
            labels = FALSE,
            k = NULL,
            split = TRUE,
            horiz = FALSE,
            ...)
  {
    # strwidth (NULL) is numeric (0), so max() returned -Inf (with a warning) on a clustering
    # without labels, and par (mar = ... + -Inf) then failed. Nothing to widen in that case.
    if (labels && (length (clustering$labels) > 0))
    {
      extra = max (graphics::strwidth (clustering$labels, units = "figure") * 30)
      if (horiz)
        opar = graphics::par (mar = graphics::par ("mar") + c (0, 0, 0, extra))
      else
        opar = graphics::par (mar = graphics::par ("mar") + c (extra, 0, 0, 0))
      on.exit (graphics::par (opar))
    }
    tree = stats::as.dendrogram (clustering)
    graphics::plot (tree, ylab = "Height",
                    leaflab = ifelse (labels, "perpendicular", "none"), horiz = horiz)
    if (is.null (k))
    {
      if (is.null (clustering$cluster))
        k = 1 + which.min (diff (sort (clustering$height, decreasing = TRUE)))
      else
        k = length (unique (clustering$cluster))
    }
    # With a single cluster there is nothing to outline -- rect.hclust() would just draw one
    # rectangle around the whole tree. This read 'split = TRUE', which looks like the
    # opposite of the intent.
    if (k == 1)
      split = FALSE
    if (split)
      # as.hclust() so that an agnes object (which the documentation accepts) works too:
      # as.dendrogram() above handles it, but rect.hclust() needs a proper hclust.
      stats::rect.hclust (stats::as.hclust (clustering), k = k, border = 2:(k + 1))
  }

#' @keywords internal
# Shared body of the print methods of the clustering classes: they all hold a 'cluster' vector
# and the data it was obtained on, and differ only by the parameters worth recalling.
printclustering <-
  function (title, x, extra = NULL)
  {
    clusters = x$cluster
    sizes = table (clusters [clusters > 0])
    noise = sum (clusters == 0)
    shortprint (title,
                c (list ("dataset" = sizelabel (x$data),
                         "clusters" = length (sizes),
                         "sizes" = paste (as.vector (sizes), collapse = ", "),
                         "unclustered" = if (noise > 0) paste (noise, "observations") else NULL),
                   extra))
  }

#' Print a DBSCAN clustering
#'
#' Prints the number of clusters found, their sizes, the observations left as noise, and the
#' two parameters used -- instead of dumping the underlying list.
#' @name print.dbs
#' @param x The clustering (object of class \code{\link{dbs-class}}, created by
#' \code{\link{DBSCAN}}).
#' @param ... Other parameters.
#' @export
#' @method print dbs
#' @seealso \code{\link{dbs-class}}, \code{\link{DBSCAN}}
#' @examples
#' require (datasets)
#' data (iris)
#' DBSCAN (iris [, -5], minpts = 5, eps = 0.65)
print.dbs <-
  function (x, ...)
    printclustering ("DBSCAN clustering", x,
                     list ("eps" = signif (x$eps, 4), "minpts" = x$MinPts))

#' Print an EM clustering
#'
#' Prints the number of clusters found, their sizes and the log-likelihood reached.
#' @name print.em
#' @param x The clustering (object of class \code{\link{em-class}}, created by \code{\link{EM}}).
#' @param ... Other parameters.
#' @export
#' @method print em
#' @seealso \code{\link{em-class}}, \code{\link{EM}}
#' @examples
#' require (datasets)
#' data (iris)
#' EM (iris [, -5], 3)
print.em <-
  function (x, ...)
  {
    # mclust records the size of the problem as n / d rather than keeping the data.
    if (is.null (x$data) && (!is.null (x$n)))
      x$data = matrix (0, nrow = x$n, ncol = x$d)
    printclustering ("EM clustering (Gaussian mixture)", x,
                     list ("model" = x$modelName,
                           "log-likelihood" = if (!is.null (x$loglik)) signif (x$loglik, 6) else NULL))
  }

#' Print a mean shift clustering
#'
#' Prints the number of clusters found, their sizes and the kernel used.
#' @name print.meanshift
#' @param x The clustering (object of class \code{\link{meanshift-class}}, created by
#' \code{\link{MEANSHIFT}}).
#' @param ... Other parameters.
#' @export
#' @method print meanshift
#' @seealso \code{\link{meanshift-class}}, \code{\link{MEANSHIFT}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' MEANSHIFT (iris [, -5])
#' }
print.meanshift <-
  function (x, ...)
    printclustering ("Mean shift clustering", x,
                     list ("kernel" = x$kernel,
                           "bandwidth" = paste (signif (x$bandwidth, 4), collapse = ", ")))

#' Print a spectral clustering
#'
#' Prints the number of clusters found, their sizes and the dimension of the spectral
#' projection.
#' @name print.spectral
#' @param x The clustering (object of class \code{\link{spectral-class}}, created by
#' \code{\link{SPECTRAL}}).
#' @param ... Other parameters.
#' @export
#' @method print spectral
#' @seealso \code{\link{spectral-class}}, \code{\link{SPECTRAL}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' SPECTRAL (iris [, -5], 3)
#' }
print.spectral <-
  function (x, ...)
    printclustering ("Spectral clustering", x,
                     list ("projection" = paste (ncol (x$proj), "dimensions")))

#' Print a self-organising map
#'
#' Prints the size of the map, how many of its units are actually used, and the dataset it was
#' fitted on.
#' @name print.som
#' @param x The map (object of class \code{\link{som-class}}, created by \code{\link{SOM}}).
#' @param ... Other parameters.
#' @export
#' @method print som
#' @seealso \code{\link{som-class}}, \code{\link{SOM}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' SOM (iris [, -5], 4, 4)
#' }
print.som <-
  function (x, ...)
  {
    grid = x$som$grid
    shortprint ("Self-organising map",
                list ("dataset" = sizelabel (x$data),
                      "map" = paste0 (grid$xdim, " x ", grid$ydim, " ", grid$topo, " grid"),
                      "units used" = paste (length (unique (x$cluster)), "of",
                                            grid$xdim * grid$ydim)))
  }

#' Clustering using \emph{K}-medoids (PAM)
#'
#' Partitions the data into \code{k} clusters around \emph{medoids} -- actual observations of
#' the dataset -- rather than around means. Being an observation, a medoid can be shown to
#' students as a representative example of its cluster, and the method tolerates outliers much
#' better than \emph{K}-means, which drags a mean towards them.
#' @name PAM
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of clusters.
#' @param criterion How \code{k} is chosen, as in \code{\link{KMEANS}}. With \code{"none"} (the
#' default) \code{k} is used as it is; with \code{"silhouette"} it is chosen between 2 and
#' \code{k} by the mean silhouette width, the criterion PAM itself optimises the closest.
#' @param graph A logical indicating whether the criterion curve is plotted.
#' @param seed A specified seed for random number generation. PAM's initialisation is
#' deterministic, so this only matters for the criterion search.
#' @param ... Other parameters, passed to \code{\link[cluster]{pam}}.
#' @return The clustering, as an object of class \code{pam} (see \code{\link[cluster]{pam}}),
#' with a \code{cluster} component holding the assignments and a \code{medoids} one holding the
#' representative observations.
#' @export
#' @seealso \code{\link{KMEANS}}, \code{\link[cluster]{pam}}, \code{\link{kmeans.getk}}
#' @examples
#' require (datasets)
#' data (iris)
#' model = PAM (iris [, -5], 3)
#' model$medoids
#' table (model$cluster, iris [, 5])
PAM <-
  function (d, k = 9, criterion = c ("none", "silhouette"), graph = FALSE, seed = NULL, ...)
  {
    if (!requireNamespace ("cluster", quietly = TRUE))
      stop ("PAM: the 'cluster' package is needed. Please install it.")
    setseed (seed)
    criterion = match.arg (criterion [1], c ("none", "silhouette"))
    kk = k
    if (criterion != "none")
    {
      if (k < 2)
        stop ("PAM: 'k' must be at least 2 when it is the largest number of clusters to ",
              "consider, but is ", k, ".")
      dis = stats::dist (d)
      sizes = 2:k
      measure = sapply (sizes, function (i)
        mean (cluster::silhouette (cluster::pam (d, i, ...)$clustering, dis) [, "sil_width"]))
      kk = sizes [which.max (measure)]
      if (graph)
      {
        graphics::plot (sizes, measure, t = "b", xlab = "Number of clusters",
                        ylab = "Mean silhouette width")
        graphics::abline (v = kk, lty = 3)
      }
    }
    res = cluster::pam (d, kk, ...)
    res$cluster = res$clustering
    # plotclus() and the rest of the package read cluster centres from '$centers'; PAM's are
    # its medoids.
    res$centers = res$medoids
    res$data = as.matrix (d)
    return (res)
  }

#' Predict function for PAM
#'
#' Returns the cluster of the closest medoid, for a new dataset.
#' @name predict.pam
#' @param object The clustering (created by \code{\link{PAM}}).
#' @param newdata A new dataset (a \code{data.frame}), with the same variables as the learning
#' dataset.
#' @param ... Other parameters.
#' @return A vector of cluster numbers.
#' @export
#' @method predict pam
#' @seealso \code{\link{PAM}}, \code{\link{predict.kmeans}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = PAM (d$train.x, 3)
#' table (predict (model, d$test.x), d$test.y)
predict.pam <-
  function (object, newdata, ...)
    apply (flexclust::dist2 (object$medoids, newdata), 2, which.min)

#' Predict function for a self-organising map
#'
#' Returns the cluster of the closest unit of the map, for a new dataset.
#' @name predict.som
#' @param object The map (created by \code{\link{SOM}}).
#' @param newdata A new dataset (a \code{data.frame}), with the same variables as the learning
#' dataset.
#' @param ... Other parameters.
#' @return A vector of cluster numbers.
#' @export
#' @method predict som
#' @seealso \code{\link{SOM}}, \code{\link{plot.som}}
#' @examples
#' \donttest{
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = SOM (d$train.x, 4, 4)
#' table (predict (model, d$test.x), d$test.y)
#' }
predict.som <-
  function (object, newdata, ...)
  {
    units = apply (flexclust::dist2 (object$som$codes [[1]], newdata), 2, which.min)
    # 'nodes' maps a unit of the map to the cluster it was assigned to.
    return (object$nodes [units])
  }

#' Predict function for hierarchical clustering
#'
#' Returns the cluster whose centre is closest, for a new dataset. A dendrogram says nothing
#' about observations it was not built on, so the rule is the usual one: the clusters of the
#' cut are summarised by their centres, and a new observation joins the nearest.
#' @name predict.hca
#' @param object The clustering (created by \code{\link{HCA}}).
#' @param newdata A new dataset (a \code{data.frame}), with the same variables as the learning
#' dataset.
#' @param k The number of clusters the dendrogram is cut into. Defaults to the cut
#' \code{\link{HCA}} already made, when it made one.
#' @param ... Other parameters.
#' @return A vector of cluster numbers.
#' @export
#' @method predict hca
#' @seealso \code{\link{HCA}}, \code{\link{predict.kmeans}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' model = HCA (d$train.x, k = 3, method = "ward")
#' table (predict (model, d$test.x), d$test.y)
predict.hca <-
  function (object, newdata, k = NULL, ...)
  {
    if (is.null (k))
      k = if (!is.null (object$cluster)) length (unique (object$cluster))
          else 1 + which.min (diff (sort (object$height, decreasing = TRUE)))
    clusters = stats::cutree (object, k)
    if (is.null (object$data))
      stop ("predict.hca: this clustering does not carry the dataset it was built on, so the ",
            "cluster centres cannot be computed. Rebuild it with HCA().")
    centres = t (sapply (sort (unique (clusters)),
                         function (i) colMeans (object$data [clusters == i, , drop = FALSE])))
    return (apply (flexclust::dist2 (centres, newdata), 2, which.min))
  }
