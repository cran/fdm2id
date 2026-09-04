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
