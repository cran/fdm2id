# Tests for text.R
#
# Deliberately restricted to functions that are simple counting/statistics (no model
# fitting, no network access): getvocab(), frequentwords(), and a small TEXTMINING() +
# NB() classification pipeline on a small text2vec::movie_review subset. GloVe-based
# functions (vectorize.words(), query.words()) are NOT covered here: this session already
# showed that GloVe fitting on a tiny corpus is fragile territory that needs hands-on
# verification (see NEWS.md), which isn't a good fit for an unattended test suite.

test_that ("getvocab() / frequentwords() run on a small local corpus (no network needed)", {
  skip_if_not_installed ("text2vec")
  skip_if_not_installed ("stopwords")
  skip_if_not_installed ("SnowballC")
  text = c ("paris is the capital of france",
           "berlin is the capital of germany",
           "france and germany are neighboring countries",
           "the capital of france is paris",
           "the capital of germany is berlin")
  vocab = getvocab (text, mincount = 2)
  expect_true (nrow (vocab) > 0)
  expect_true (all (c ("term", "term_count") %in% colnames (vocab)))
  top = frequentwords (text, 3, mincount = 1)
  expect_true (length (top) <= 3)
})

test_that ("TEXTMINING + NB classification pipeline runs on a small movie_review subset", {
  skip_if_not_installed ("text2vec")
  skip_if_not_installed ("e1071")
  data ("movie_review", package = "text2vec")
  d = movie_review [1:60, 2:3]
  d [, 1] = factor (d [, 1])
  d = splitdata (d, 1, seed = 0)
  model = TEXTMINING (d$train.x, NB, labels = d$train.y, mincount = 5)
  pred = predict (model, d$test.x)
  # d$test.x is a plain character vector here (the lone "review" column collapses via
  # [s, -target] indexing), which is exactly what TEXTMINING()'s corpus argument expects.
  expect_equal (length (pred), length (d$test.x))
})

# --- Regression: TEXTMINING (vector = "words") produced an unusable object ---------------
# The "words" branch built a "textmining" object with no vectorizer slot, while
# predict.textmining() dereferences it unconditionally.

test_that ("predict() explains why a word-based textmining model cannot predict documents", {
  skip_if_not_installed ("text2vec")
  skip_if_not_installed ("stopwords")
  set.seed (1)
  vocabulary = paste0 ("w", 1:200)
  corpus = replicate (200, paste (sample (vocabulary, 40, replace = TRUE), collapse = " "))
  model = suppressMessages (TEXTMINING (corpus, HCA, vector = "words",
                                        k = 3, mincount = 5, ndim = 5, maxiter = 3))
  expect_equal (model$vector, "words")
  expect_null (model$vectorizer)
  expect_error (predict (model, corpus [1:5]), "vector = \"words\"")
})

# --- Regression: vectorize.docs() always densified the document-term matrix ---------------

test_that ("vectorize.docs (sparse = TRUE) returns the same values, sparsely", {
  skip_if_not_installed ("text2vec")
  skip_if_not_installed ("stopwords")
  set.seed (1)
  vocabulary = paste0 ("w", 1:300)
  corpus = replicate (300, paste (sample (vocabulary, 50, replace = TRUE), collapse = " "))
  dense = vectorize.docs (corpus = corpus, mincount = 2, transform = "tfidf")
  sparse = vectorize.docs (corpus = corpus, mincount = 2, transform = "tfidf", sparse = TRUE)
  expect_s3_class (dense, "data.frame")
  expect_false (is.data.frame (sparse))
  expect_equal (dim (dense), dim (sparse))
  expect_equal (unname (as.matrix (dense)), unname (as.matrix (sparse)))
  expect_lt (as.numeric (utils::object.size (sparse)), as.numeric (utils::object.size (dense)))
})

# =========================================================================================
# Seventh batch of the audit: CRAN policy and hygiene
# =========================================================================================

# --- Regression: plotzipf() left options(scipen=) set in the user's session -----------------

test_that ("plotzipf() restores the options it changes", {
  skip_if_not_installed ("text2vec")
  set.seed (1)
  words = c ("science", "data", "model", "cluster", "random",
             "forest", "tree", "vector", "matrix", "kernel")
  corpus = replicate (20, paste (sample (words, 40, replace = TRUE, prob = 1 / (1:10)),
                                 collapse = " "))
  before = getOption ("scipen")
  grDevices::pdf (NULL)
  on.exit (grDevices::dev.off ())
  plotzipf (corpus)
  expect_identical (getOption ("scipen"), before)
})

# --- Regression: loadtext() downloaded into the user's home directory -----------------------

test_that ("loadtext() defaults to the session's temporary directory", {
  expect_identical (deparse (formals (loadtext)$dir), "tempdir()")
})

test_that ("loadtext() reads a plain file and a zip archive from an explicit directory", {
  f = tempfile (fileext = ".txt")
  writeLines (c ("alpha beta", "gamma delta"), f)
  expect_equal (loadtext (f), "alpha beta gamma delta")
  expect_equal (loadtext (f, collapse = FALSE), c ("alpha beta", "gamma delta"))
  skip_if (Sys.which ("zip") == "")
  src = file.path (tempdir (), "loadtext-src")
  dir.create (src, showWarnings = FALSE)
  writeLines ("one two", file.path (src, "a.txt"))
  writeLines ("three four", file.path (src, "b.txt"))
  archive = file.path (tempdir (), "loadtext-corpus.zip")
  unlink (archive)
  old = setwd (src)
  on.exit (setwd (old))
  utils::zip (archive, c ("a.txt", "b.txt"), flags = "-q")
  setwd (old)
  # A directory that does not exist yet is created rather than failing.
  target = file.path (tempdir (), "loadtext-out", "nested")
  unlink (target, recursive = TRUE)
  expect_equal (sort (loadtext (archive, dir = target)), c ("one two", "three four"))
  expect_true (dir.exists (target))
})
