# Tests for factorial.R
# FactoMineR is a Depends of fdm2id (always installed alongside it), so no
# skip_if_not_installed() guard is needed here.

test_that ("PCA smoke test + kaiser rule", {
  data (iris)
  pca = PCA (iris [, 1:4])
  expect_s3_class (pca, "factorial")
  k = kaiser (pca)
  expect_true (k >= 1 && k <= 4)
})

test_that ("PCA with a qualitative supplementary variable", {
  data (iris)
  pca = PCA (iris, quali.sup = 5)
  expect_s3_class (pca, "factorial")
  expect_output (print (pca))
})

test_that ("CA / MCA smoke tests", {
  data (children, package = "FactoMineR")
  ca = CA (children, row.sup = 15:18, col.sup = 6:8)
  expect_s3_class (ca, "factorial")
  data (tea, package = "FactoMineR")
  mca = MCA (tea, quanti.sup = 19, quali.sup = 20:36)
  expect_s3_class (mca, "factorial")
})

test_that ("plot.factorial runs for individual/eigenvalue PCA plots without error", {
  data (iris)
  pca = PCA (iris, quali.sup = 5)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (plot (pca), NA)
  expect_error (plot (pca, type = "eig"), NA)
})

# =========================================================================================
# Sixth batch of the audit
# =========================================================================================

# --- Regression: plot.factorial() dropped 'col' and 'labels' for CA and MCA -----------------

test_that ("plot.factorial() forwards col and labels for CA and MCA too", {
  skip_if_not_installed ("FactoMineR")
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  d = discretizeDF (iris [, -5],
                    default = list (method = "interval", breaks = 3,
                                    labels = c ("small", "medium", "large")))
  mca = MCA (d)
  expect_error (plot (mca, col = as.numeric (iris [, 5]), labels = FALSE), NA)
  expect_error (plot (mca, labels = TRUE), NA)
  ca = CA (as.data.frame.matrix (table (d [, 1], d [, 2])))
  expect_error (plot (ca, col = c (1, 2, 3), labels = TRUE), NA)
})

# --- Regression: the 'factorial' S4 declaration is gone, the objects are unchanged -----------

test_that ("factorial results still answer is() and carry their sub-class", {
  data (iris)
  p = PCA (iris [, -5])
  expect_true (methods::is (p, "factorial"))
  expect_true (inherits (p, "pca"))
})

# =========================================================================================
# Tenth batch: projecting new observations into a factorial space
# =========================================================================================

test_that ("predict.factorial() gives FactoMineR's own supplementary-individual coordinates", {
  skip_if_not_installed ("FactoMineR")
  data (iris)
  d = splitdata (iris, 5, seed = 0)
  pca = PCA (d$train.x)
  both = rbind (d$train.x, d$test.x)
  ref = FactoMineR::PCA (both, ind.sup = nrow (d$train.x) + seq_len (nrow (d$test.x)),
                         ncp = ncol (both), graph = FALSE)
  expect_equal (unname (predict (pca, d$test.x)), unname (ref$ind.sup$coord))
  # The axes are those of the training analysis; the new observations do not move them.
  expect_equal (unname (pca$ind$coord), unname (ref$ind$coord))
})

test_that ("an observation of the training set projects onto its own coordinates", {
  data (iris)
  pca = PCA (iris [, -5])
  expect_equal (unname (predict (pca, iris [1:5, -5])), unname (pca$ind$coord [1:5, ]),
                tolerance = 1e-10)
})

test_that ("predict.factorial() works for MCA and CA too", {
  skip_if_not_installed ("FactoMineR")
  data (iris)
  d = discretizeDF (iris [, -5],
                    default = list (method = "interval", breaks = 3,
                                    labels = c ("s", "m", "l")))
  mca = MCA (d [1:100, ])
  # The reference is asked for the same number of axes: fdm2id keeps every one the analysis can
  # produce, where FactoMineR's own default stops at five.
  ref = FactoMineR::MCA (d, ind.sup = 101:150, ncp = ncol (mca$ind$coord), graph = FALSE)
  expect_equal (unname (predict (mca, d [101:150, ])), unname (ref$ind.sup$coord))
  tab = as.data.frame.matrix (table (d [, 1], d [, 2]))
  extra = rbind (tab, extra = c (3, 5, 2))
  ca = CA (tab)
  refc = FactoMineR::CA (extra, row.sup = nrow (extra), ncp = ncol (ca$row$coord), graph = FALSE)
  expect_equal (unname (predict (ca, extra [nrow (extra), , drop = FALSE])),
                unname (refc$row.sup$coord))
})

test_that ("predict.factorial() does not need the supplementary variables, and checks names", {
  data (iris)
  pca = PCA (iris, quali.sup = 5)
  # 'Species' is supplementary: it plays no part in the axes, so it need not be supplied.
  expect_equal (nrow (predict (pca, iris [1:5, -5])), 5)
  expect_equal (unname (predict (pca, iris [1:5, -5])), unname (pca$ind$coord [1:5, ]),
                tolerance = 1e-10)
  expect_error (predict (pca, data.frame (Unknown = 1:3)), "does not know")
  expect_error (predict (structure (list (), class = "factorial"), iris [1:3, -5]),
                "not built by fdm2id")
})

# --- Fifth audit: the Kaiser rule counts, it does not index -------------------------------

test_that ("kaiser() answers a number of axes, never -Inf", {
  data (iris)
  expect_equal (kaiser (PCA (iris, quali.sup = 5)), 1)
  # max (which (...)) over an empty set is max (integer (0)), i.e. -Inf with a warning: that
  # is what uncorrelated standardised variables produce, every axis carrying the average.
  fake = list (eig = cbind (rep (1, 4), rep (25, 4), cumsum (rep (25, 4))))
  expect_equal (kaiser (fake), 0)
  expect_true (kaiser (PCA (iris [, -5])) >= 1)
})

# =========================================================================================
# Twelfth audit: ncp defaults to every axis the analysis can produce
# =========================================================================================

test_that ("MCA keeps every axis it can produce", {
  data (zoo)
  data (credit)
  data (titanic)
  # A variable with J_q modalities carries J_q - 1 axes, so Q active variables give J - Q.
  # The course material used to reach it by passing ncp = 100 by hand.
  for (case in list (list (zoo, list (quali.sup = 17)), list (credit, list ()),
                     list (titanic, list ())))
  {
    model = do.call (MCA, c (list (case [[1]]), case [[2]]))
    expect_equal (ncol (model$ind$coord), nrow (model$eig))
  }
  # The zoo practical asks how many axes an MCA produces, and answers twenty.
  expect_equal (nrow (MCA (zoo, quali.sup = 17)$eig), 20)
  expect_equal (ncol (MCA (zoo, quali.sup = 17)$ind$coord), 20)
  # An explicit ncp is still honoured, in both directions.
  expect_equal (ncol (MCA (credit, ncp = 3)$ind$coord), 3)
  expect_equal (ncol (MCA (zoo, quali.sup = 17, ncp = 100)$ind$coord), 20)
})

test_that ("the modalities are counted as observed, not as declared", {
  data (titanic)
  # On the first 200 rows the declared levels promise six axes and the analysis produces two:
  # counting the declared ones would ask FactoMineR for axes that do not exist.
  expect_equal (fdm2id:::mca.ncp (titanic [1:200, ]), 2)
  expect_equal (nrow (MCA (titanic [1:200, ])$eig), 2)
  expect_equal (fdm2id:::mca.ncp (titanic), 6)
  # Supplementary individuals and variables take no part in the count.
  expect_equal (fdm2id:::mca.ncp (titanic, quali.sup = 4),
                sum (sapply (titanic [, -4], nlevels) - 1))
})

test_that ("CA keeps every axis it can produce", {
  data (children, package = "FactoMineR")
  data (titanic)
  model = CA (children, row.sup = 15:18, col.sup = 6:8)
  expect_equal (ncol (model$row$coord), nrow (model$eig))
  expect_equal (nrow (model$eig), fdm2id:::ca.ncp (children, row.sup = 15:18, col.sup = 6:8))
  # min (I, J) - 1, the first dimension of a contingency table being its margins.
  table = as.data.frame.matrix (table (titanic [, 1], titanic [, 4]))
  expect_equal (fdm2id:::ca.ncp (table), 1)
  expect_equal (nrow (CA (table)$eig), 1)
})

test_that ("PCA counts the observations as well as the variables", {
  data (iris)
  data (cookies.desc.train)
  expect_equal (fdm2id:::pca.ncp (iris, quali.sup = 5), 4)
  # 40 observations of 700 variables: the number of axes is bounded by n - 1, not by p.
  expect_equal (fdm2id:::pca.ncp (cookies.desc.train), 39)
  expect_equal (ncol (PCA (cookies.desc.train)$ind$coord), 39)
  expect_equal (fdm2id:::pca.ncp (iris [, -5], ind.sup = 1:100), 4)
})

test_that ("an analysis with no axis to compute says so", {
  data (titanic)
  # Every variable takes a single value on those twenty rows; FactoMineR fails inside its
  # singular value decomposition on "max(nu, nv) must be positive".
  expect_error (MCA (titanic [1:20, ]), "no axis to compute")
  expect_error (CA (data.frame (a = c (5, 3, 2))), "no axis to compute")
})
