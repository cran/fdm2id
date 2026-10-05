# fdm2id 1.0.2

## Residuals

* `resplot (index = "fitted")` plots the studentized residuals against the fitted values, the
  plot that was missing: no value of `index` gave it. `index = 0`, which plots them against the
  response, now warns -- a least-squares residual is correlated with the response
  (`sqrt (1 - R^2)`) whatever the model, so that plot shows a trend that says nothing about the
  fit. The examples used it; they now use `"fitted"`. New `xlab` argument, and a vector given as
  `index` names the axis.
* `residuals()`, `fitted()` and `rstudent()` apply to a model built by the package, without
  going through `model$model`.

## Bugs

* `LR()` learnt from a one-column data frame can now predict: the variable kept its own name
  when learning, and was renamed `X` when predicting.
* `treeplot()`: the rectangles take the colour of their cluster in `plotclus()` (its
  `cutree()` number + 1); they were coloured from left to right of the tree, so the
  dendrogram and the scatter plot of a clustering did not match.
* `evaluation.kappa()` returned `NA` beyond 46 340 predictions (`n * n` overflowed as an
  integer), as repeated bootstraps easily reach.
* `plotclus()` and `scatterplot()` placed their legend by the first two columns of the data
  rather than by the plotted coordinates -- the PCA plane, beyond two variables -- and so often
  on the points.
* `DBSCAN (epsilonDist)`, `EM (clusters)` and `params` in `SVR()`, `SVRl()`, `SVRr()` and
  `MLPREG()`, renamed in 1.0.0, were swallowed by `...` without a word: an old script ran with
  the defaults and gave another result. They are now an error that names the new argument.

## Graphics

* Legends are placed where they hide the least of the plot: of the eight positions `legend()`
  offers along the edges, the one covering the fewest points -- or the fewest stretches of
  curve, bar or box -- and, at equal count, the one furthest from them. `legendpos = "auto"`
  is the new default of `plotdata()`, `plot()` on a factorial analysis, a CDA or a feature
  selection, and `boxclus()`; `scatterplot (legend = "auto")` replaces `"auto1"` and
  `"auto2"`, still accepted. An explicit position is honoured. The legends that had a fixed
  place -- ROC and cost curves, `kmeans.getk (graph = TRUE)`, the PCR/PLS and ridge/lasso
  tuning curves, `plotzipf()`, the scree plot -- are placed the same way.
* New `plot()` method for `HCA()`: the dendrogram, the branches of each cluster in the colour
  `plotclus()` gives it, those above the cut dashed.
* `treeplot (highlight = "branches")` colours the branches of the dendrogram below the cut
  with the colour of their cluster -- the colours of `plotclus()`, so that the dendrogram and
  the scatter plot of the same clustering match -- and dashes the ones above it, all in the
  current line width (`par (lwd = 2)` for thick lines). `ylab` is now an argument, and the
  other arguments are passed to `plot()` (`axes = FALSE`).
* `plot()` on a `SOM()`, `type = "mapping"`: the cells are coloured with a quarter of the
  colour of their cluster instead of two thirds, so that the points, in the full colour of the
  same cluster, stand out.
* The legends drawn by the package follow `par ("cex.lab")`, as the axis labels do: unchanged
  by default, and enlarged with them (`par (cex.lab = 1.5)`).

## Changes

* `selectfeatures()` and `FEATURESELECTION()`: the default multivariate criterion is CFS
  rather than mRMR. With `algorithm = "ranking"`, the criterion scores the nested subsets of
  the best-ranked features; mRMR, `fstat` and `inertiaratio` are averages, which only decrease
  as less relevant features join, so the default always kept a single feature.
* `kmeans.getk()` re-seeds before each number of clusters, so the partition behind each value
  is the one `KMEANS (d, k, nstart, seed)` returns, and the retained `k` the one the pseudo-F of
  those partitions gives.
* `plotdata()`: new argument `scale` (`FALSE` by default, as before), to project scaled
  variables with `type = "pca"` and `"scatter"`.

## Documentation

* `LINREG (quali)`: the four values described.
* `performance (methodparameters)`: the tuning done beforehand when it is missing, its cost,
  and how to pass pre-tuned parameters. `BAGGING()` and `ADABOOST()`: how to give the base
  learner through `performance()`.

## Packaging

* `datasets` is no longer imported in the `NAMESPACE`. The package code never used it -- only
  the examples, the tests and one vignette do, and they attach it themselves -- and R-devel
  now notes a base package that is imported while listed under `Suggests`.

# fdm2id 1.0.1

* The vignettes are pre-computed: their code is run when the sources are prepared, by
  `make-vignettes.R`, rather than every time the package is built. Rebuilding them took
  between six and eleven minutes on the CRAN check machines -- most of the check time, and
  more than its budget. It now costs a pandoc pass. The results they print are unchanged.
* Shorter introductions in the vignettes.

# fdm2id 1.0.0

The package was reviewed with Claude Code. Some forty bugs were fixed as a result, and a test
suite was added.

## Breaking changes

Results that change: `evaluation.precision()`, `evaluation.recall()`, `evaluation.kappa()`,
`evaluation.adjr2()`, `compare.jaccard()`, `intern.dunn()`, `CDA()`, `MLPREG()`, `LINREG()`,
`TSNE()`, `correlated()`, `predict()` on `SVM`, `performance (type = "roc" / "cost")`, feature
selection (CFS, Relief, mRMR), and every split or cross-validation, now stratified
(`stratify`). `performance()` also re-seeds before each method it is given, so a method that
searches a grid no longer scores differently depending on its position in the list.

Arguments:

* the final arguments of every learning method are now `nfolds`, `tune`, `methodparameters`,
  `graph`, `seed`, in that order
* `HCA()`: second argument
* `DBSCAN (epsilonDist)` is now `eps`; `EM (clusters)` is now `k`; `params` is now
  `methodparameters` in `SVR()`, `SVRl()`, `SVRr()` and `MLPREG()`
* new defaults: `graph = FALSE` everywhere, `STUMP (randomvar)`,
  `FEATURESELECTION (multieval)`, `LINREG (regeval, nrep, validation)`, `CART (xval)`,
  `loadtext (dir)`, `performance (fuzzy)`, `CA (ncp)`, `MCA (ncp)`

Return values: `compare (comp = "cluster")`, `LINREG (reg = "subset")`, `MLPREG (tune = TRUE)`,
`FEATURESELECTION (tune = TRUE)`.

Errors raised where a value was returned: misspelt `comp`, `type`, `quali` and `reg`, `CDA()`,
`LR()`, `intern.dunn()`, `performance()`.

Removed: `HDBSCAN()`, `ISOFOREST()`, `UMAP()`, the `setClass()` declarations.

Datasets: `ozone` has 13 variables instead of 10; `movies` is now simulated and rated from
1 to 10.

## New features

* new methods: `GBREG()`, `PAM()`, penalized logistic regression (`LR (reg = )`)
* `predict()` on factorial analyses and on every clustering method
* multiclass precision, recall and F-measure (`average`)
* `performance()` on a test set given directly, without resampling
* number of clusters: silhouette, gap statistic, elbow
* `selectfeatures (unieval = "randomforest")` and `plot()` on a selection
* `plotdata (type = "correlation")` and `plotdata (target = )`
* `HCA (engine = )`, automatic `eps` in `DBSCAN()` and number of clusters in `HCA()`
* `seed` in every learning method, `k` in every clustering method
* package objects print a summary
* `legendpos` in `plot()` on a `CDA` object
* faster hyperparameter search, `SPECTRAL()`, `HCA()` and feature selection
* eight vignettes, one per case study of the courses, listed by
  `vignette (package = "fdm2id")`; they print the same figures as the course handouts, and
  flag in comments every call whose result moves from one run to the next without a `seed`
