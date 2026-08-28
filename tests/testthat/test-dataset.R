# Tests for dataset.R (synthetic dataset generators)
# All generators default to graph = TRUE; passing graph = FALSE avoids plotdata() entirely,
# so no null graphics device is needed here.

test_that ("synthetic dataset generators produce data.frames of the expected shape", {
  d1 = data.gauss (n = 50, graph = FALSE, seed = 0)
  expect_equal (nrow (d1), 50)
  expect_true ("Class" %in% colnames (d1))

  d2 = data.parabol (graph = FALSE, seed = 0)
  expect_true (nrow (d2) > 0)
  expect_true ("Class" %in% colnames (d2))

  d3 = data.diag (n = 50, graph = FALSE, seed = 0)
  expect_equal (nrow (d3), 50)

  d4 = data.target1 (graph = FALSE, seed = 0)
  expect_true (nrow (d4) > 0)

  d5 = data.target2 (graph = FALSE, seed = 0)
  expect_true (nrow (d5) > 0)

  d6 = data.twomoons (graph = FALSE, seed = 0)
  expect_true (nrow (d6) > 0)

  d7 = data.xor (n = 20, graph = FALSE, seed = 0)
  expect_true (nrow (d7) > 0)
  expect_true ("Class" %in% colnames (d7))
})

test_that ("synthetic dataset generators are reproducible given the same seed", {
  d1 = data.gauss (n = 50, graph = FALSE, seed = 0)
  d2 = data.gauss (n = 50, graph = FALSE, seed = 0)
  expect_equal (d1, d2)

  d3 = data.twomoons (graph = FALSE, seed = 0)
  d4 = data.twomoons (graph = FALSE, seed = 0)
  expect_equal (d3, d4)
})

# --- Regression: the data.* generators honour a seed set by the caller ------------------
# They all used to call set.seed (NULL), which re-seeds from the clock (see test-misc.R).

test_that ("the data.* generators are reproducible under a caller-set seed", {
  generators = list ("data.gauss"    = function () data.gauss (n = 50, graph = FALSE),
                     "data.twomoons" = function () data.twomoons (n = 30, graph = FALSE),
                     "data.parabol"  = function () data.parabol (n = c (30, 10), graph = FALSE),
                     "data.target1"  = function () data.target1 (n = 20, graph = FALSE),
                     "data.diag"     = function () data.diag (n = 30, graph = FALSE),
                     "data.xor"      = function () data.xor (n = 20, graph = FALSE))
  for (name in names (generators))
  {
    set.seed (42)
    a = generators [[name]] ()
    set.seed (42)
    b = generators [[name]] ()
    expect_equal (a, b, info = name)
  }
})

# --- Regression: data.target1()'s sigma was documented but ignored ------------------------

test_that ("data.target1() honours sigma", {
  a = data.target1 (n = 400, sigma = 0.01, graph = FALSE, seed = 1)
  b = data.target1 (n = 400, sigma = 0.50, graph = FALSE, seed = 1)
  expect_false (isTRUE (all.equal (a [, 1], b [, 1])))
  # The radial spread must grow with sigma.
  spread = function (d) stats::sd (sqrt (d [, 1]^2 + d [, 2]^2) -
                                     as.numeric (as.character (factor (d [, 3], labels = 1:3))))
  expect_lt (spread (a), spread (b))
})

# --- Regression: data.xor() had a duplicated assignment to 'levels' -------------------------

test_that ("data.xor() honours custom class labels", {
  d = data.xor (n = 10, graph = FALSE, seed = 0)
  expect_equal (levels (d$Class), c ("Class 1", "Class 2"))
  d = data.xor (n = 10, graph = FALSE, seed = 0, levels = c ("neg", "pos"))
  expect_equal (levels (d$Class), c ("neg", "pos"))
})

# --- The movies dataset is simulated, and stays that way ----------------------------------

test_that ("movies is a simulated 1-10 ratings matrix", {
  data (movies)
  expect_true (is.matrix (movies))
  expect_equal (dim (movies), c (49, 55))
  expect_true (is.integer (movies))
  expect_true (all (movies >= 1 & movies <= 10))
  expect_false (anyNA (movies))
  # The film titles are real, the viewers are not.
  expect_true ("Star Wars (1977)" %in% rownames (movies))
  expect_true (all (grepl ("^User ", colnames (movies))))
  # A 1-5 scale would mean the MovieLens ratings had come back.
  expect_gt (max (movies), 5)
})

test_that ("the structure the visualisation methods read is preserved", {
  data (movies)
  # The films everyone rates highly are still the ones that were rated highly, and the first
  # principal axis still carries about a quarter of the variance, as in the source data.
  variance = 100 * stats::prcomp (movies)$sdev [1]^2 / sum (stats::prcomp (movies)$sdev^2)
  expect_gt (variance, 20)
  expect_lt (variance, 35)
  # Star Wars keeps the company it kept.
  d = as.matrix (stats::dist (movies))
  diag (d) = Inf
  expect_true ("Empire Strikes Back, The (1980)" %in%
               rownames (movies) [order (d ["Star Wars (1977)", ]) [1:3]])
})
