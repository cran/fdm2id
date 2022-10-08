#' DBSCAN model
#'
#' This class contains the model obtained by the DBSCAN method.
#' @name dbs-class
#' @slot cluster A vector of integers indicating the cluster to which each point is allocated.
#' @slot eps Reachability distance (parameter).
#' @slot MinPts Reachability minimum no. of points (parameter).
#' @slot isseed A logical vector indicating whether a point is a seed (not border, not noise).
#' @slot data The dataset that has been used to fit the map (as a \code{matrix}).
#' @exportClass dbs
#' @seealso \code{\link{DBSCAN}}
setClass ("dbs",
          representation (cluster = "vector",
                          eps = "numeric",
                          MinPts = "numeric",
                          isseed = "vector",
                          data = "matrix"))

#' Expectation-Maximization model
#'
#' This class contains the model obtained by the EM method.
#' @name em-class
#' @slot modelName A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
#' @slot prior Specification of a conjugate prior on the means and variances.
#' @slot n The number of observations in the dataset.
#' @slot d The number of variables in the dataset.
#' @slot G The number of components of the mixture.
#' @slot z A matrix whose \code{[i,k]}th entry is the conditional probability of the ith observation belonging to the kth component of the mixture.
#' @slot parameters A names list giving the parameters of the model.
#' @slot control A list of control parameters for EM.
#' @slot loglik The log likelihood for the data in the mixture model.
#' @slot cluster A vector of integers (from \code{1:k}) indicating the cluster to which each point is allocated.
#' @exportClass em
#' @seealso \code{\link{EM}}, \code{\link[mclust]{mclustModelNames}}
setClass ("em",
          representation (modelName = "character",
                          prior = "numeric",
                          n = "numeric",
                          d = "numeric",
                          G = "numeric",
                          z = "matrix",
                          parameters = "list",
                          control = "list",
                          loglik = "numeric",
                          cluster = "vector"))

#' MeanShift model
#'
#' This class contains the model obtained by the MEANSHIFT method.
#' @name meanshift-class
#' @slot cluster A vector of integers indicating the cluster to which each point is allocated.
#' @slot value A vector or matrix containing the location of the classified local maxima in the support.
#' @slot data The leaning set.
#' @slot kernel A string indicating the kernel associated with the kernel density estimate that the mean shift is optimizing over.
#' @slot bandwidth Used in the kernel density estimate for steepest ascent classification.
#' @slot alpha A scalar tuning parameter for normal kernels.
#' @slot iterations The number of iterations to perform mean shift.
#' @slot epsilon A scalar used to determine when to terminate the iteration of a individual query point.
#' @slot epsilonCluster A scalar used to determine the minimum distance between distinct clusters.
#' @exportClass meanshift
#' @seealso \code{\link{MEANSHIFT}}
setClass ("meanshift",
          representation (cluster = "vector",
                          value = "vector",
                          data = "matrix",
                          kernel = "character",
                          bandwidth = "vector",
                          alpha = "numeric",
                          iterations = "numeric",
                          epsilon = "numeric",
                          epsilonCluster = "numeric"))

#' Self-Organizing Maps model
#'
#' This class contains the model obtained by the SOM method.
#' @name som-class
#' @slot som An object of class \code{kohonen} representing the fitted map.
#' @slot nodes A \code{vector} of integer indicating the cluster to which each node is allocated.
#' @slot cluster A \code{vector} of integer indicating the cluster to which each observation is allocated.
#' @slot data The dataset that has been used to fit the map (as a \code{matrix}).
#' @exportClass som
#' @seealso \code{\link{plot.som}}, \code{\link{SOM}}, \code{\link[kohonen]{som}}
setClass ("som",
          representation (som = "list",
                          nodes = "vector",
                          cluster = "vector",
                          data = "matrix"))

#' Spectral clustering model
#'
#' This class contains the model obtained by Spectral clustering.
#' @name spectral-class
#' @slot cluster A \code{vector} of integer indicating the cluster to which each observation is allocated.
#' @slot proj The projection of the dataset in the spectral space.
#' @slot centers The cluster centers (on the spectral space).
#' @exportClass spectral
#' @seealso \code{\link{SPECTRAL}}
setClass ("spectral",
          representation (cluster = "vector",
                          proj = "matrix",
                          centers = "matrix"))

#' @keywords internal
accuracy0 <-
  function (clus, gt)
  {
    return (sum (diag (table (clus, gt))) / length (clus))
  }

#' @keywords internal
accuracy1 <-
  function (clus, gt)
  {
    if (sum (clus) == 0)
      return (NA)
    res = accuracy0 (clus, gt == min (gt))
    for (i in (min (gt) + 1):(max (gt)))
    {
      tmp = accuracy0 (clus, gt == i)
      if (tmp > res)
        res = tmp
    }
    return (res)
  }

#' Clustering Box Plots
#'
#' Produce a box-and-whisker plot for clustering results.
#' @name boxclus
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param clusters Cluster labels of the training set (\code{vector} or \code{factor}).
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
    noise = min (clusters) == 0
    mini = min (d)
    maxi = max (d)
    nbclusters = length (unique (clusters))
    names = colnames (d)
    nbclusters * ncol (d)
    at = (0:(ncol (d) - 1)) * nbclusters + 1 + (nbclusters - 1) / 2
    v = (1:(ncol (d) - 1)) * nbclusters + .5
    d = utils::stack (d)
    d$cluster = clusters
    graphics::boxplot (values~cluster+ind, d, ylim = c (mini, maxi), ylab = "", col = 2:(nbclusters + 1), xaxt='n', xlab = "")
    graphics::axis (side = 1, at = at, labels = names, lwd.ticks = FALSE, lwd = 0)
    graphics::abline (v = v, lty = 2, col = "grey")
    labels = sort (unique (clusters))
    if (noise)
      labels = c ("Noise", paste ("Cluster", labels [-1]))
    else
      labels = paste ("Cluster", labels)
    graphics::legend (x = legendpos, legend = labels, fill = sort (unique (1 + clusters)), bty = "n")
  }

#' Comparison of two sets of clusters
#'
#' Comparison of two sets of clusters
#' @name compare
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @param eval The evluation criterion.
#' @param comp Indicates whether a "max" or a "pairwise" evaluation should be used, or the evaluation for each individual "cluster".
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
    res = NULL
    for (e in eval)
    {
      tmp = get (paste ("compare", e, sep = ".")) (clus, gt, comp)
      res = c (res, tmp)
    }
    if (!is.null (res))
      names (res) = eval
    return (res)
  }

#' Comparison of two sets of clusters, using accuracy
#'
#' Comparison of two sets of clusters, using accuracy
#' @name compare.accuracy
#' @param clus The extracted clusters.
#' @param gt The real clusters.
#' @param comp Indicates whether a "max" or a "pairwise" evaluation should be used, or the evaluation for each individual "cluster".
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
    kk1 = as.numeric (clus)
    kk2 = as.numeric (gt)
    res = 0
    if (comp [1] == "pairwise")
    {
      comparisonmatrix = comparison.matrix (length (kk1))
      p1 = comparisonmatrix %*% kk1 != 0
      p2 = comparisonmatrix %*% kk2 != 0
      t = table (p1, p2)
      res = sum (diag (t)) / sum (t)
    }
    else
    {
      res = NULL
      clusters = min (kk1):max (kk1)
      for (i in clusters)
        res = c (res, accuracy1 (kk1 == i, kk2))
      if (comp [1] == "max")
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
#' @param comp Indicates whether a "max" or a "pairwise" evaluation should be used, or the evaluation for each individual "cluster".
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
    kk1 = as.numeric (clus)
    kk2 = as.numeric (gt)
    res = 0
    if (comp [1] == "pairwise")
    {
      comparisonmatrix = comparison.matrix (length (kk1))
      p1 = comparisonmatrix %*% kk1 != 0
      p2 = comparisonmatrix %*% kk2 != 0
      t = table (p1, p2)
      res = sum (diag (t)) / sum (t)
    }
    else
    {
      res = NULL
      clusters = min (kk1):max (kk1)
      for (i in clusters)
        res = c (res, jaccard1 (kk1 == i, kk2))
      if (comp [1] == "max")
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
#' @param comp Indicates whether a "max" or a "pairwise" evaluation should be used, or the evaluation for each individual "cluster".
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
    kk1 = as.numeric (clus)
    kk2 = as.numeric (gt)
    res = 0
    if (comp [1] == "pairwise")
    {
      comparisonmatrix = comparison.matrix (length (kk1))
      p1 = comparisonmatrix %*% kk1 != 0
      p2 = comparisonmatrix %*% kk2 != 0
      res = irr::kappa2 (cbind (p1, p2), weight = "equal")$value
    }
    else
    {
      res = NULL
      clusters = min (kk1):max (kk1)
      for (i in clusters)
        res = c (res, kappa1 (kk1 == i, kk2))
      if (comp [1] == "max")
        res = stats::weighted.mean (res, table (kk1), na.rm = TRUE)
      else
        names (res) = paste ("Cluster", clusters)
    }
    return (res)
  }

#' @keywords internal
comparison.matrix <-
  function (n)
  {
    res = NULL
    if (n >= 2)
    {
      for (i in (n - 1):1)
      {
        a = matrix (0, nrow = i, ncol = i)
        diag (a) = -1
        if (n - i - 1 > 0)
          res = rbind (res, cbind (matrix(0, nrow = i, ncol = (n - i - 1)), rep (1, i), a))
        else
          res = rbind (res, cbind (rep (1, i), a))
      }
    }
    else
      stop ("n should be bigger than 2")
    return (res)
  }

#' DBSCAN clustering method
#'
#' Run the DBSCAN algorithm for clustering.
#' @name DBSCAN
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param minpts Reachability minimum no. of points.
#' @param epsilonDist Reachability distance.
#' @param ... Other parameters.
#' @return A clustering model obtained by DBSCAN.
#' @export
#' @seealso \code{\link[fpc]{dbscan}}, \code{\link{dbs-class}}, \code{\link{distplot}}, \code{\link{predict.dbs}}
#' @examples
#' require (datasets)
#' data (iris)
#' DBSCAN (iris [, -5], minpts = 5, epsilonDist = 1)
DBSCAN <-
  function (d, minpts, epsilonDist, ...)
  {
    res = fpc::dbscan (d, MinPts = minpts, eps = epsilonDist)
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
    Kdistance = sort (apply (as.matrix (stats::dist (d, upper = TRUE, diag = TRUE)), 2, sort) [k + 1, ], decreasing = TRUE)
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
#' @param clusters Either an integer (the number of clusters) or a (\code{vector}) indicating the cluster to which each point is initially allocated.
#' @param model A character string indicating the model. The help file for \code{\link[mclust]{mclustModelNames}} describes the available models.
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
  function (d, clusters, model = "VVV", ...)
  {
    if (length (clusters) == 1)
      clusters = stats::kmeans (d, clusters, nstart = 10)$cluster
    z = mclust::unmap (clusters)
    p = mclust::mstep (data = d, modelName = model, z = z)
    res = mclust::em (data = d, modelName = model, parameters = p$parameters)
    res = c (res, list (cluster = apply (res$z, 1, which.max)))
    class (res) = "em"
    return (res)
  }

#' Hierarchical Cluster Analysis method
#'
#' Run the HCA method for clustering.
#' @name HCA
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param method Character string defining the clustering method.
#' @param k The number of cluster.
#' @param ... Other parameters.
#' @return The cluster hierarchy (\code{hca} object).
#' @export
#' @seealso \code{\link[cluster]{agnes}}
#' @examples
#' require (datasets)
#' data (iris)
#' HCA (iris [, -5], method = "ward", k = 3)
HCA <-
  function (d, method = c ("ward", "single"), k = NULL, ...)
  {
    hc = cluster::agnes (d, method = method [1])
    if (is.null (k))
      k = 1 + which.min (diff (sort (hc$height, decreasing = TRUE)))
    cluster = stats::cutree (hc, k)
    hc = stats::as.hclust (hc)
    r = c (hc, list (cluster = cluster, k = k))
    class (r) = c ("hca", class (hc))
    return (r)
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
#' intern (km$clus, iris [, -5])
#' intern (km$clus, iris [, -5], type = "cluster")
#' intern (km$clus, iris [, -5], eval = c ("intraclass", "interclass"))
#' intern (km$clus, iris [, -5], eval = c ("intraclass", "interclass"), type = "cluster")
intern <-
  function (clus, d, eval = "intraclass", type = c ("global", "cluster"))
  {
    res = sapply (eval, function (e)
    {
      return (get (paste ("intern", e, sep = ".")) (clus, d, type))
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
#' @param type Indicates whether a "global" or a "cluster"-wise evaluation should be used.
#' @return The evaluation of the clustering.
#' @export
#' @seealso \code{\link{intern}}, \code{\link{intern.interclass}}, \code{\link{intern.intraclass}}
#' @examples
#' require (datasets)
#' data (iris)
#' km = KMEANS (iris [, -5], k = 3)
#' intern.dunn (km$clus, iris [, -5])
intern.dunn <-
  function (clus, d, type = c ("global"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    if (type [1] != "global")
    {
      message ("Dunn index only works for global evaluation.")
      return (NULL)
    }
    dis = flexclust::dist2 (d, d)
    clusters = sort (unique (clus))
    dmin = min (sapply (clusters, function (i) sapply (clusters, function (j)
    {
      if (i != j)
        return (min (dis [clus == i, clus == j]))
      else
        return (NA)
    })), na.rm = TRUE)
    dmax = max (sapply (clusters, function (i)
    {
      return (max (dis [clus == i, clus == i]))
    }))
    return (dmin / dmax)
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
#' intern.interclass (km$clus, iris [, -5])
intern.interclass <-
  function (clus, d, type = c ("global", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
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
#' intern.intraclass (km$clus, iris [, -5])
intern.intraclass <-
  function (clus, d, type = c ("global", "cluster"))
  {
    if (!is.vector (clus))
      clus = clus$cluster
    centers = apply (d, 2, function (v) tapply (v, clus, mean))
    indices = sort (unique (clus))
    res = sapply (indices, function (index)
    {
      center = matrix (centers [index, ], nrow = 1)
      return (sum (flexclust::dist2 (center, d [clus == index, ])^2))
    })
    if (type [1] == "global")
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
    res = jaccard0 (clus, gt == min (gt))
    for (i in (min (gt) + 1):(max (gt)))
    {
      tmp = jaccard0 (clus, gt == i)
      if (tmp > res)
        res = tmp
    }
    return (res)
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
    res = kappa0 (clus, gt == min (gt))
    for (i in (min (gt) + 1):(max (gt)))
    {
      tmp = kappa0 (clus, gt == i)
      if (tmp > res)
        res = tmp
    }
    return (res)
  }

#' K-means method
#'
#' Run K-means for clustering.
#' @name KMEANS
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of cluster.
#' @param criterion The criterion for cluster number selection. If \code{none}, \code{k} is used, if not the number of cluster is selected between 2 and \code{k}.
#' @param graph A logical indicating whether or not a graphic should be plotted (cluster number selection).
#' @param nstart Define how many random sets should be chosen.
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
  function (d, k = 9, criterion = c ("none", "pseudo-F"), graph = FALSE, nstart = 10, ...)
  {
    kk = k
    if (criterion [1] == "pseudo-F")
      kk = kmeans.getk (d = d, max = k, criterion = criterion, graph = graph, nstart = nstart)
    return (stats::kmeans (d, kk, nstart = nstart))
  }

#' Estimation of the number of clusters for \emph{K}-means
#'
#' Estimate the optimal number of cluster of the \emph{K}-means clustering method.
#' @name kmeans.getk
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param max The maximum number of clusters. Values from 2 to \code{max} are evaluated.
#' @param criterion The criterion to be optimized. \code{"pseudo-F"} is the only criterion implemented in the current version.
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param nstart The number of random sets chosen for \code{\link[stats]{kmeans}} initialization.
#' @param seed A specified seed for random number generation.
#' @return The optimal number of cluster of the \emph{K}-means clustering method according to the chosen criterion.
#' @export
#' @seealso \code{\link{pseudoF}}, \code{\link[stats]{kmeans}}
#' @examples
#' require (datasets)
#' data (iris)
#' kmeans.getk (iris [, -5])
kmeans.getk <-
  function (d, max = 9, criterion = "pseudo-F", graph = TRUE, nstart = 10, seed = NULL)
  {
    set.seed (seed)
    k = NA
    measure = vector ("numeric", max - 1)
    measure2 = NULL
    criterion2 = NULL
    if (criterion == "pseudo-F")
    {
      measure2 = vector ("numeric", max - 1)
      criterion2 = as.expression (substitute (R^2))
      for (i in 2:max)
      {
        km = stats::kmeans (d, i, nstart = nstart)
        measure [i - 1] = pseudoF (km)
        measure2 [i - 1] = km$betweenss / km$totss
      }
      k = which.max (measure) + 1
    }
    if (graph & !is.na (k))
    {
      if (!is.null (measure2))
      {
        opar = graphics::par (mar = c(5, 4, 4, 5) + .1)
        on.exit (graphics::par (opar))
      }
      graphics::plot (2:max, measure, t = "b", xlab = "Number of clusters", ylab = criterion)
      if (!is.null (measure2))
      {
        opar = graphics::par (new = TRUE)
        on.exit (graphics::par (opar))
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
#' @param epsilon A scalar used to determine when to terminate the iteration of a individual query point.
#' @param epsilonCluster A scalar used to determine the minimum distance between distinct clusters.
#' @param ... Other parameters.
#' @return The clustering (\code{meanshift} object).
#' @export
#' @seealso \code{\link[meanShiftR]{meanShift}}, \code{\link{predict.meanshift}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' MEANSHIFT (iris [, -5], bandwidth = .75)
#' }
MEANSHIFT <-
  function (d, mskernel = "NORMAL", bandwidth = rep (1, ncol (d)), alpha = 0, iterations = 10, epsilon = 1e-08, epsilonCluster = 1e-04, ...)
  {
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
      centers.coord = x$som$codes
      xlab = colnames (d) [1]
      ylab = colnames (d) [2]
      col = x$cluster + 1
      if (ncol (x$data) > 2)
      {
        pca = FactoMineR::PCA (d, scale.unit = FALSE, graph = FALSE, ncp = 2)
        centers.coord = t ((t (x$som$codes [[1]]) - apply (d, 2, mean)) / apply (d, 2, stats::sd))
        centers.coord = centers.coord %*% pca$svd$V
        dd = pca$ind$coord
        row.names (dd) = row.names (d)
        d = dd
        xlab = paste ("Prin. 1 (", round (pca$eig [1, 2], 2), " %)", sep = "")
        ylab = paste ("Prin. 2 (", round (pca$eig [2, 2], 2), " %)", sep = "")
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
      adj [adj [, 1] < adj [, 2], ]
      graphics::segments (centers.coord [adj [, 1], 1], centers.coord [adj [, 1], 2],
                          centers.coord [adj [, 2], 1], centers.coord [adj [, 2], 2],
                          lty = 2, col = "darkgrey")
      adj = which (as.matrix (stats::dist (pts [, 1])) ==  0 & as.matrix (stats::dist (pts [, 2])) ==  1, arr.ind = TRUE)
      adj [adj [, 1] < adj [, 2], ]
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
        bgcol [empty] = col [apply (flexclust::dist2(d, x$som$codes [[1]] [empty, ]), 2, which.min)]
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
      graphics::plot (x$som, type = "mapping", classif = x$unit.classif,
                      pch = 19, col = pcol, bgcol = bgcol, main = "", labels = lab)
    }
    else
      message (paste (type, ": Unknown"))
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
#' ward = HCA (iris [, -5], method = "ward", k = 3)
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
    method = class (clustering)
    clusters = NULL
    centres = NULL
    if ("kmeans" %in% method)
    {
      clusters = clustering$cluster
      if (centers)
        centres = clustering$centers
    }
    else if ("fclust" %in% method)
    {
      clusters = clustering$cluster
      if (centers)
        centres = clustering$centers
    }
    else if ("em" %in% method)
    {
      clusters = apply (clustering$z, 1, which.max)
      if (centers)
        centres = t (clustering$parameters$mean)
    }
    else if ("dbs" %in% method)
    {
      clusters = clustering$cluster
    }
    else if ("som" %in% method)
    {
      clusters = clustering$cluster
    }
    else if ("meanshift" %in% method)
    {
      clusters = clustering$cluster
    }
    else if ("spectral" %in% method)
    {
      clusters = clustering$cluster
    }
    else if ("hca" %in% method)
    {
      if (is.null (k))
      {
        if (is.null (clustering$cluster))
          k = 1 + which.min (diff (sort (clustering$height, decreasing = TRUE)))
        else
          k = length (unique (clustering$cluster))
      }
      clusters = stats::cutree (clustering, k)
    }
    if ((type [1] == "tree") & ("hca" %in% method))
      treeplot (clustering, k = k, ...)
    else if (type [1] == "height" & ("hca" %in% method))
      graphics::barplot (utils::tail (sort (clustering$height), n = tailsize), names.arg = tailsize:1,
                         xlab = "Number of clusters", ylab = "Height", main = "", sub = "")
    else if ((type [1] == "scatter") & (!is.vector (d)) & ("som" %in% method))
      plot.som (clustering, type = type, ...)
    else if ((type [1] == "mapping") & ("som" %in% method))
      plot.som (clustering, type = type, ...)
    else if (type [1] == "scatter")
      scatterplot (d, clusters, centres, ...)
    else if (type [1] == "boxplot")
      boxclus (d, clusters, ...)
    else if (type [1] == "words")
      plotcloud (d, clusters, ...)
    else
      stop ("Uncorrect parameters. Plot type is not available for the clustering method.")
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
#' @param clusters Cluster labels of the training set (\code{vector} or \code{factor}).
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
    noise = min (clusters) == 0
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
      pca = FactoMineR::PCA (d, scale.unit = FALSE, graph = FALSE, ncp = 2)
      dd = pca$ind$coord
      row.names (dd) = row.names (d)
      if (!is.null (centers))
      {
        m = apply (centers, 2, mean)
        centers = sweep (centers, 2, m)
        centers = centers %*% pca$svd$V
      }
      xlab = paste ("Prin. 1 (", round (pca$eig [1, 2], 2), " %)", sep = "")
      ylab = paste ("Prin. 2 (", round (pca$eig [2, 2], 2), " %)", sep = "")
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
      z = mclust::unmap (clusters [clusters != 0])
      p = mclust::mstep (data = dd [clusters != 0, ], modelName = "VVV", z = z)
      if (is.null (centers))
        centers = t (p$parameters$mean)
      for (i in unique (clusters))
        if (i != 0)
        {
          col = i + 1
          e = car::ellipse (p$parameters$variance$sigma [,, i], centre = centers [i, ], npoints = 1000)
          graphics::polygon (e, border = col, lty = 2)
          eig = eigen (p$parameters$variance$sigma [,, i])
          seg = sweep (eig$vectors, 1, sqrt (eig$values), FUN = "*")
          seg = seg * sqrt (stats::qchisq (0.95, 2))
          graphics::segments (seg [1, 1] + centers [i, 1], -seg [1, 2] + centers [i, 2],
                              -seg [1, 1] + centers [i, 1], seg [1, 2] + centers [i, 2],
                              col = col, lty = 3)
          graphics::segments (seg [2, 1] + centers [i, 1], -seg [2, 2] + centers [i, 2],
                              -seg [2, 1] + centers [i, 1], seg [2, 2] + centers [i, 2],
                              col = col, lty = 3)
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
    labels = sort (unique (clusters))
    if (noise)
      labels = c ("Noise", paste ("Cluster", labels [-1]))
    else
      labels = paste ("Cluster", labels)
    graphics::legend (x = legend, legend = labels, col = sort (unique (col)), pch = 1)
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
#' @return The fitted Kohonen's map as an object of class \code{som}.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{plot.som}}, \code{\link{som-class}}, \code{\link[kohonen]{som}}
#' @examples
#' require (datasets)
#' data (iris)
#' SOM (iris [, -5], xdim = 5, ydim = 5, post = "ward", k = 3)
SOM <-
  function (d, xdim = floor (sqrt (nrow (d))), ydim = floor (sqrt (nrow (d))), rlen = 10000, post = c ("none", "single", "ward"), k = NULL, ...)
  {
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

#' Spectral clustering method
#'
#' Run a Spectral clustering algorithm.
#' @name SPECTRAL
#' @param d The dataset (\code{matrix} or \code{data.frame}).
#' @param k The number of cluster.
#' @param sigma Width of the gaussian used to build the affinity matrix.
#' @param graph A logical indicating whether or not a graphic should be plotted (projection on the spectral space of the affinity matrix).
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{spectral-class}}
#' @examples
#' \dontrun{
#' require (datasets)
#' data (iris)
#' SPECTRAL (iris [, -5], k = 3)
#' }
SPECTRAL <-
  function (d, k, sigma = 1, graph = TRUE, ...)
  {
    a = exp (-flexclust::dist2 (d, d)^2 / (2 * sigma * sigma))
    diag (a) = 0
    p = diag (1 / sqrt (rowSums (a)))
    l = p %*% a %*% p
    x = eigen (l)$vectors [, 1:k]
    proj = sweep (x, 1, sqrt (rowSums (x * x)), "/")
    colnames (proj) = paste ("Comp.", 1:k)
    rownames (proj) = rownames (d)
    km = stats::kmeans (proj, centers = k, nstart = 100)
    cluster = km$cluster
    if (graph)
      plotdata (d = proj, k = factor (paste ("Cluster", cluster)))
    res = list (cluster = cluster, proj = proj, centers = km$centers)
    class (res) = "spectral"
    return (res)
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
        methodNames = as.character (match.call ()$clusteringmethods)
        if (length (methodNames) > 1)
          methodNames = methodNames [-1]
        if (length (methodNames) != length (clusteringmethods))
          methodNames = NULL
      }
    }
    clusteringmethods = c (clusteringmethods)
    if (is.null (originals))
    {
      if (!is.null (seed))
        set.seed (seed + 1)
      else
        set.seed (seed)
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
        set.seed (seed)
        s = matrix (sample (nrow (d), nsampling * nrow (d), replace = T), ncol = nsampling)
        rescomp = apply (s, 2, function (v)
        {
          clusters = clus (d [v, ], graph = graph, ...)
          if (graph)
            plotclus (clusters, d [v, ])
          compare (original [unique (v)], clusters$cluster [!duplicated (v)], eval = eval, comp = comp)
        })
        rescomp = apply (s, 2, function (v) compare (original [unique (v)], clus (d [v, ], graph = graph,...)$cluster [!duplicated (v)], eval = eval, comp = comp))
        if (is.vector (rescomp))
        {
          res = mean (rescomp, na.rm = TRUE)
          names (res) = eval
        }
        else if ((is.matrix (rescomp)) & (length (eval) == 1))
        {
          res = matrix (rescomp, ncol = nsampling)
          res = matrix (apply (res, 1, mean, na.rm = TRUE), ncol = 1)
          colnames (res) = eval
          rownames (res) = clusternames
        }
        else if (type [1] == "cluster")
        {
          rescomp = unlist (rescomp)
          dim = c (length (rescomp) / (nsampling * length (eval)), length (eval), nsampling)
          res = array (rescomp, dim)
          res = apply (res, 1:2, mean, na.rm = TRUE)
          colnames (res) = eval
          rownames (res) = clusternames
        }
        else
        {
          res = apply (rescomp, 1, mean, na.rm = TRUE)
          names (res) = eval
        }
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
#' hca = HCA (iris [, -5], method = "ward", k = 3)
#' treeplot (hca)
treeplot <-
  function (clustering,
            labels = FALSE,
            k = NULL,
            split = TRUE,
            horiz = FALSE,
            ...)
  {
    if (labels)
    {
      if (horiz)
      {
        opar = graphics::par (mar = graphics::par ("mar") + c (0, 0, 0, max (graphics::strwidth (clustering$labels, units = "figure") * 30)))
        on.exit (graphics::par (opar))
      }
      else
      {
        opar = graphics::par (mar = graphics::par("mar") + c (max (graphics::strwidth (clustering$labels, units = "figure") * 30), 0, 0, 0))
        on.exit (graphics::par (opar))
      }
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
    if (split)
      stats::rect.hclust (clustering, k = k, border = 2:(k + 1))
  }
