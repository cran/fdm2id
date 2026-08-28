# Tests for visualization.R

test_that ("plotdata() runs without error for the common plot types", {
  data (iris)
  grDevices::pdf (NULL) # discard graphics output, avoid writing Rplots.pdf during tests
  on.exit (grDevices::dev.off ())
  expect_error (plotdata (iris [, -5], iris [, 5]), NA)
  expect_error (plotdata (iris [, -5], iris [, 5], type = "scatter"), NA)
  expect_error (plotdata (iris [, -5], iris [, 5], type = "boxplot"), NA)
})

test_that ("SVD returns ind/var projections of the expected shape", {
  data (iris)
  res = SVD (iris [, -5])
  expect_equal (nrow (res$proj$ind), nrow (iris))
  expect_equal (nrow (res$proj$var), ncol (iris [, -5]))
})

test_that ("TSNE is reproducible given a seed", {
  skip_if_not_installed ("Rtsne")
  data (iris)
  res1 = TSNE (iris [, -5], seed = 0, nstart = 1)
  res2 = TSNE (iris [, -5], seed = 0, nstart = 1)
  expect_equal (res1$Y, res2$Y)
})

# --- Regression: plotdata() silently drew nothing for an unsupported type ----------------
# plotdata.matrix() had no final else, unlike plotdata.vector().

test_that ("plotdata() reports an unsupported type on a multi-variable dataset", {
  data (iris)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  expect_message (plotdata (iris [, -5], iris [, 5], type = "words"), "unavailable")
  expect_message (plotdata (iris [, -5], iris [, 5], type = "not_a_type"), "unavailable")
  # A supported type must of course stay silent.
  expect_silent (plotdata (iris [, -5], iris [, 5], type = "boxplot"))
})

# ========================================================================================
# Regression tests for the "batch 4" fixes.
# ========================================================================================

# --- Regression: TSNE()'s nstart loop and its de-duplicated output -----------------------
# The loop compared the current best with its own cost -- always FALSE -- so the nstart - 1
# extra embeddings were computed and thrown away. And Rtsne needs distinct rows, so the
# result had fewer rows than x whenever the data contained a repeated observation (iris has
# one), which silently mismatched points and colours in plotdata (type = "tsne").

test_that ("TSNE() returns one row per observation, duplicates included", {
  skip_if_not_installed ("Rtsne")
  data (iris)
  expect_lt (nrow (unique (iris [, -5])), nrow (iris)) # iris does contain a duplicate
  res = TSNE (iris [, -5], perplexity = 10, nstart = 1, seed = 1)
  expect_equal (nrow (res$Y), nrow (iris))
  expect_equal (length (res$costs), nrow (iris))
  # The two identical observations must land on the same point.
  keys = do.call (paste, c (as.data.frame (iris [, -5]), sep = "\r"))
  i = which (duplicated (keys)) [1]
  j = which (keys == keys [i]) [1]
  expect_equal (unname (res$Y [i, ]), unname (res$Y [j, ]))
})

test_that ("TSNE() actually keeps the best of its nstart runs", {
  skip_if_not_installed ("Rtsne")
  data (iris)
  set.seed (7)
  one = utils::tail (TSNE (iris [, -5], perplexity = 10, nstart = 1)$itercosts, 1)
  set.seed (7)
  many = utils::tail (TSNE (iris [, -5], perplexity = 10, nstart = 5)$itercosts, 1)
  expect_lte (many, one)
})

# --- Regression: plotdata()'s alpha and asp were accepted but ignored ---------------------

test_that ("plotdata() honours alpha and asp", {
  data (iris)
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  expect_error (plotdata (iris [, -5], iris [, 5], type = "scatter", alpha = 100), NA)
  expect_error (plotdata (iris [, -5], iris [, 5], type = "scatter", alpha = 255), NA)
  expect_error (plotdata (iris [, -5], iris [, 5], type = "scatter", asp = NA), NA)
  # addalpha() is vectorised and produces 8-digit hex colours.
  cols = fdm2id:::addalpha (c (1, 2, 3), 200)
  expect_equal (length (cols), 3)
  expect_true (all (grepl ("^#[0-9A-Fa-f]{8}$", cols)))
})

# --- Regression: histogram()'s density overlay was clipped to ylim = c (0, 1) -------------

test_that ("histogram() does not clip a density above 1", {
  grDevices::pdf (file = tempfile (fileext = ".pdf"))
  on.exit (grDevices::dev.off ())
  set.seed (1)
  x = stats::rnorm (1000, 0, 0.1) # peaks around 4
  expect_gt (max (stats::density (x)$y), 1)
  expect_error (fdm2id:::histogram (x, ""), NA)
})


# =========================================================================================
# Third audit, E1 / E2
# =========================================================================================

test_that ("SVD projects on a single axis", {
  data (iris)
  # diag() of a *single* number builds an identity matrix of that size, so ndim = 1 failed
  # instead of projecting on the first axis.
  res = SVD (iris [, -5], ndim = 1)
  expect_equal (dim (res$proj$ind), c (nrow (iris), 1))
  expect_equal (dim (res$proj$var), c (4, 1))
  # And it is the right projection, not merely a non-error: rank-one reconstruction.
  expect_equal (as.matrix (res$proj$ind) %*% t (res$v),
                as.matrix (iris [, -5]) %*% (res$v %*% t (res$v)),
                ignore_attr = TRUE)
})

test_that ("SVD names the dimensions it cannot deliver", {
  data (iris)
  expect_error (SVD (iris [, -5], ndim = 0), "between 1 and 4")
  expect_error (SVD (iris [, -5], ndim = 99), "between 1 and 4")
  expect_error (SVD (iris [, -5], ndim = 2.5), "whole number")
})

test_that ("plotdata (type = 'cda') works on a two-class target", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  two = iris
  levels (two [, 5]) = c ("+", "+", "-")
  # Two classes give a single canonical axis, and proj [, 1:2] asked for a column that does
  # not exist ("undefined columns selected").
  expect_error (plotdata (two [, -5], two [, 5], type = "cda"), NA)
  expect_error (plotdata (two [, -5], two [, 5], type = "cda", labels = TRUE), NA)
  expect_error (plotdata (iris [, -5], iris [, 5], type = "cda"), NA)
  # The single axis is named, as plot.cda() has always assumed.
  expect_equal (colnames (CDA (two [, -5], two [, 5])$proj), "Can. 1")
})

# --- Ninth audit: too many variables for a matrix of scatter plots -------------------------

test_that ("plotdata (type = 'pairs') says what to do when the panels do not fit", {
  data (movies)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # 55 variables means 55 x 55 panels: R stopped on "figure margins too large", which says
  # nothing about the alternatives.
  expect_error (plotdata (movies, type = "pairs"), "do not fit on this device")
  expect_error (plotdata (movies, type = "pairs"), "type = \"pca\"")
  # And a dataset that does fit is unaffected.
  data (iris)
  expect_error (plotdata (iris, type = "pairs"), NA)
})

# --- Ninth audit: '...' goes to the method or to the plot, not to both ---------------------

test_that ("plotdata routes an argument to the method that understands it", {
  data (alcohol)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # 'perplexity' is TSNE()'s, not graphics::plot()'s, which used to answer
  # '"perplexity" is not a graphical parameter' once per element it drew. The tsne branch
  # papered over it with suppressWarnings(); the arguments are now routed instead.
  expect_silent (plotdata (alcohol, type = "tsne", labels = TRUE, perplexity = 5))
  expect_silent (plotdata (alcohol, type = "tsne", perplexity = 5))
  # A graphical argument still reaches the plot, alongside a method one.
  expect_silent (plotdata (alcohol, type = "tsne", perplexity = 5, main = "title"))
  # And the method really receives it: a perplexity too large for the sample is refused by
  # Rtsne, which it could not be if the argument were being swallowed on the way.
  expect_error (plotdata (alcohol, type = "tsne", perplexity = 30), "perplexity")
  # Same routing for the map.
  expect_silent (plotdata (alcohol, type = "som", labels = TRUE, xdim = 4, ydim = 4))
})

test_that ("the routing keeps method and plot arguments apart", {
  route = fdm2id:::plotdata.route (list (perplexity = 5, main = "t", nstart = 2, col = 3), TSNE)
  expect_named (route$method, c ("perplexity", "nstart"))
  expect_named (route$plot, c ("main", "col"))
  # Nothing named goes nowhere, and nothing goes to both.
  expect_equal (length (route$method) + length (route$plot), 4)
  expect_equal (fdm2id:::plotdata.route (list (), TSNE), list (method = list (), plot = list ()))
})

# =========================================================================================
# Eleventh audit: plotdata (type = "correlation")
# =========================================================================================

test_that ("a numeric target gives the signed correlation, strongest first", {
  data (mtcars)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  res = plotdata (mtcars [, -1], mtcars [, 1], type = "correlation")
  expect_named (res)
  expect_setequal (names (res), colnames (mtcars) [-1])
  # The values are Pearson's, sign included.
  for (v in names (res))
    expect_equal (unname (res [v]), unname (stats::cor (mtcars [[v]], mtcars [, 1])), info = v)
  # Sorted by strength: barplot draws the first bar at the bottom, so the weakest comes first.
  expect_equal (abs (res), sort (abs (res)))
  expect_equal (names (res) [length (res)], "wt")
})

test_that ("a categorical target gives the correlation ratio", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # More than two classes: a correlation is not defined, the correlation ratio is.
  res = plotdata (iris, type = "correlation")   # Species found on its own
  expect_setequal (names (res), colnames (iris) [-5])
  expect_true (all (res >= 0 & res <= 1))
  expect_equal (unname (res ["Petal.Length"]),
                unname (sqrt (fdm2id:::fseval.inertiaratio (iris [, -5], iris [, 5]) ["Petal.Length"])))
  # Exactly two classes: the ratio is |r| with the classes coded 0/1, so the sign comes back
  # and says which class the variable is larger in.
  two = iris
  levels (two [, 5]) = c ("setosa", "other", "other")
  signed = plotdata (two [, -5], two [, 5], type = "correlation")
  expect_lt (signed ["Sepal.Width"], 0)     # setosa has the wider sepals
  expect_gt (signed ["Petal.Length"], 0)
  expect_equal (abs (signed ["Sepal.Width"]),
                sqrt (fdm2id:::fseval.inertiaratio (two [, -5], two [, 5]) ["Sepal.Width"]))
})

test_that ("the correlation plot copes with what it cannot correlate", {
  data (mtcars)
  data (titanic)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # No target to relate the variables to.
  expect_message (plotdata (mtcars, type = "correlation"), "needs numeric variables")
  # A constant variable correlates with nothing rather than giving NA.
  flat = mtcars
  flat$flat = 1
  res = plotdata (flat [, -1], flat [, 1], type = "correlation")
  expect_equal (unname (res ["flat"]), 0)
  expect_false (anyNA (res))
  # No numeric variable at all.
  expect_message (plotdata (titanic [, -4], titanic [, 4], type = "correlation"))
  # A single predictor still draws.
  expect_error (plotdata (mtcars [, "wt", drop = FALSE], mtcars [, 1], type = "correlation"), NA)
  # And the device is left as it was found, despite the widened margin.
  before = graphics::par (c ("mar", "mfrow"))
  plotdata (mtcars [, -1], mtcars [, 1], type = "correlation")
  expect_equal (graphics::par (c ("mar", "mfrow")), before)
})

test_that ("a continuous target is not turned into as many classes as observations", {
  data (mtcars)
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  # 'k' is coerced to a factor wherever it colours or groups; type = "correlation" reads it as
  # the target instead, and a numeric one would otherwise become a 32-level factor.
  expect_error (plotdata (mtcars [, -1], mtcars [, 1], type = "correlation"), NA)
  # The types that do colour by it keep working on numeric class codes -- cluster numbers are
  # the usual case, and they must stay legible.
  km = KMEANS (iris [, -5], k = 3, seed = 0)
  expect_error (plotdata (iris [, -5], km$cluster, type = "scatter"), NA)
  expect_error (plotdata (iris [, -5], km$cluster, type = "boxplot"), NA)
  expect_error (plotdata (iris [, -5], km$cluster, type = "cda"), NA)
  # And a categorical target still reaches the correlation plot.
  expect_error (plotdata (iris, type = "correlation"), NA)
})

# --- plotdata (target): the variable to explain, given apart from the colouring one --------
# 'k' has to be categorical, so a numeric variable could not be used as the target of
# type = "correlation" without being turned into as many levels as there are observations.

test_that ("plotdata (target) accepts a numeric target and gives Pearson's correlation", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  res = plotdata (iris [, 1:3], target = iris [, 4], type = "correlation")
  ref = sapply (iris [, 1:3], function (v) stats::cor (v, iris [, 4]))
  expect_equal (unname (res [names (ref)]), unname (ref))
  expect_true (any (res < 0))               # le signe est conserve
})

test_that ("plotdata (target) and plotdata (k) agree on a categorical target", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_equal (plotdata (iris [, -5], target = iris [, 5], type = "correlation"),
                plotdata (iris [, -5], iris [, 5], type = "correlation"))
})

test_that ("the correlation ratio does not depend on the order of the classes", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  shuffled = factor (iris [, 5], levels = rev (levels (iris [, 5])))
  expect_equal (plotdata (iris [, -5], target = shuffled, type = "correlation"),
                plotdata (iris [, -5], target = iris [, 5], type = "correlation"))
  eta = plotdata (iris [, -5], target = iris [, 5], type = "correlation")
  expect_true (all (eta >= 0) && all (eta <= 1))   # un rapport de correlation vit dans [0, 1]
})

test_that ("a categorical target also colours the plot, a numeric one does not", {
  data (iris)
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  expect_error (plotdata (iris [, -5], target = iris [, 5], type = "pairs"), NA)
  expect_error (plotdata (iris [, 1:3], target = iris [, 4], type = "pairs"), NA)
})
