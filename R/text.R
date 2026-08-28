#' Document vectorization object
#'
#' This class contains a vectorization model for textual documents.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{vectorizer}}{The vectorizer.}
#'   \item{\code{transform}}{The transformation to be applied after vectorization (normalization, TF-IDF).}
#'   \item{\code{phrases}}{The phrase detection method.}
#'   \item{\code{tfidf}}{The TF-IDF transformation.}
#'   \item{\code{lsa}}{The LSA transformation.}
#'   \item{\code{tokens}}{The token from the original document.}
#' }
#' @name vectorizer-class
#' @seealso \code{\link{vectorize.docs}}, \code{\link{query.docs}}
NULL

#' Text mining object
#'
#' Object used for text mining.
#'
#' Objects of this class are plain lists with the following components:
#' \describe{
#'   \item{\code{vectorizer}}{The vectorizer.}
#'   \item{\code{vectors}}{The vectorized dataset.}
#'   \item{\code{res}}{The result of the text mining method.}
#' }
#' @name textmining-class
#' @seealso \code{\link{TEXTMINING}}, \code{\link{vectorize.docs}}
NULL

#' @keywords internal
addphrases <-
  function (it, mincount = 50, maxiter = 10)
  {
    vocab = text2vec::create_vocabulary (it, stopwords = stopwords::stopwords ("en"))
    vocab = text2vec::prune_vocabulary (vocab, term_count_min = mincount)
    model = text2vec::Collocations$new(vocabulary = vocab, collocation_count_min = mincount, pmi_min = 0)
    model$fit (it)
    nphrases = 0
    iter = 0
    while ((nphrases != nrow (model$collocation_stat)) && (iter < maxiter))
    {
      iter = iter + 1
      nphrases = nrow (model$collocation_stat)
      model$prune (pmi_min = 8, gensim_min = 10, lfmd_min = -25)
      model$partial_fit (it)
    }
    return (model)
  }

#' @keywords internal
cleanup <-
  function (corpus, removesinglechars = TRUE)
  {
    res = sapply (corpus, function (text) tolower (text))
    res = sapply (res, function (text) gsub ("[^[:alnum:]]", " ", text))
    if (removesinglechars)
      res = sapply (res, function (text) gsub ("\\b[[:alnum:]]{1}\\b", "", text))
    res = sapply (res, function (text) gsub ("\\s+", " ", text))
    return (res)
  }

#' @keywords internal
createiterator <-
  function (corpus, lang,  minphrasecount = NULL, removesinglechars = TRUE)
  {
    it = tokens (corpus, lang = lang, removesinglechars = removesinglechars)
    phrases = NULL
    if ((!is.null (minphrasecount)) && (minphrasecount > 0))
    {
      phrases = addphrases (it, mincount = minphrasecount)
      it = phrases$transform (it)
    }
    return (it)
  }

#' @keywords internal
createvectorizer <-
  function (corpus, it = NULL, phrases = NULL, vocab = NULL, lang, stopwords = lang, excludewords = NULL, ngram = 1, mincount = 10, minphrasecount = NULL,
            transform = c ("none", "l1", "tfidf", "lsa"), latentdim = 50, removesinglechars = TRUE)
  {
    if (is.null (it))
      it = createiterator (corpus, lang, minphrasecount, removesinglechars = removesinglechars)
    if (is.null (vocab))
      vocab = getvocab (corpus, mincount, minphrasecount, ngram, stopwords, excludewords = excludewords, it = it, lang = lang)
    vectorizer = text2vec::vocab_vectorizer (vocab)
    res = list (vectorizer = vectorizer, transform = transform [1], minphrasecount = minphrasecount, tokens = it, phrases = phrases)
    if (transform [1] == "tfidf")
    {
      dtm = text2vec::create_dtm (it, vectorizer)
      tfidf = text2vec::TfIdf$new()
      dtm = text2vec::fit_transform(dtm, tfidf)
      res$tfidf = tfidf
    }
    else if (transform [1] == "lsa")
    {
      dtm = text2vec::create_dtm (it, vectorizer)
      tfidf = text2vec::TfIdf$new()
      dtm = text2vec::fit_transform(dtm, tfidf)
      res$tfidf = tfidf
      lsa = text2vec::LSA$new(n_topics = latentdim)
      dtm = text2vec::fit_transform(dtm, lsa)
      res$lsa = lsa
    }
    class (res) = "vectorizer"
    return (res)
  }

#' Frequent words
#'
#' Most frequent words of the corpus.
#' @name frequentwords
#' @param corpus The corpus of documents (a vector of characters) or the vocabulary of the documents (result of function \code{getvocab}).
#' @param nb The number of words to be returned.
#' @inheritParams getvocab
#' @return The most frequent words of the corpus.
#' @export
#' @seealso \code{\link{getvocab}}
#' @examples
#' data (capitals)
#' frequentwords (capitals, 10, mincount = 2)
#' vocab = getvocab (capitals, mincount = 2)
#' frequentwords (vocab, 10)
frequentwords <-
  function (corpus, nb, mincount = 5, minphrasecount = NULL, ngram = 1, lang = "en", stopwords = lang, excludewords = NULL, removesinglechars = TRUE)
  {
    vocab = NULL
    if ("text2vec_vocabulary" %in% class (corpus))
      vocab = corpus
    else
      vocab = getvocab (corpus, mincount = mincount, minphrasecount = minphrasecount, ngram = ngram, lang = lang, stopwords = stopwords, excludewords = excludewords, removesinglechars = removesinglechars)
    return (vocab [vocab [, "term_count"] >= vocab [nrow (vocab) + 1 - nb, "term_count"], "term"])
  }

#' Extract words and phrases from a corpus
#'
#' Extract words and phrases from a corpus of documents.
#' @name getvocab
#' @param corpus The corpus of documents (a vector of characters).
#' @param mincount Minimum word count to be considered as frequent.
#' @param minphrasecount Minimum collocation of words count to be considered as frequent.
#' @param ngram maximum size of n-grams.
#' @param lang The language of the documents (NULL if no stemming).
#' @param stopwords The language whose stop words are removed (\code{"en"}, ...), or
#' \code{NULL} to keep them. A list of words of your own goes to \code{excludewords}.
#' @param excludewords An optional custom vector of additional words to exclude from the vocabulary (e.g. corpus-specific stop words), on top of (or instead of) the language stopwords given through \code{stopwords}.
#' @param removesinglechars Whether single-character tokens are removed during cleanup.
#' @param ... Other parameters.
#' @return The vocabulary used in the corpus of documents.
#' @export
#' @seealso \code{\link{plotzipf}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{create_vocabulary}}
#' @examples
#' data (capitals)
#' vocab1 = getvocab (capitals, mincount = 2) # With stemming
#' nrow (vocab1)
#' vocab2 = getvocab (capitals, mincount = 2, lang = NULL) # Without stemming
#' nrow (vocab2)
#' # Excluding additional, corpus-specific words
#' vocab3 = getvocab (capitals, mincount = 2, excludewords = c ("capital", "europe"))
getvocab <-
  function (corpus, mincount = 5, minphrasecount = NULL, ngram = 1, lang = "en", stopwords = lang, excludewords = NULL, removesinglechars = TRUE, ...)
  {
    dots = list (...)
    it = NULL
    if (!is.null (dots$it))
      it = dots$it
    else
    {
      it = tokens (corpus, lang = lang, removesinglechars = removesinglechars)
      if ((!is.null (minphrasecount)) && (minphrasecount > 0))
      {
        phrases = addphrases (it, mincount = minphrasecount)
        it = phrases$transform (it)
      }
    }
    sw = character(0)
    if (!is.null (stopwords))
      sw = stopwords::stopwords (stopwords)
    if (!is.null (excludewords))
      sw = union (sw, excludewords)
    vocab = text2vec::create_vocabulary (it, ngram = c (1, ngram), stopwords = sw)
    vocab = text2vec::prune_vocabulary (vocab, term_count_min = mincount)
    return (vocab)
  }

#' load a text file
#'
#' (Down)Load a text file (and extract it if it is in a zip file).
#' @name loadtext
#' @param file The path or URL of the text file. If not specified, defaults to an interactive file chooser (\code{\link[base]{file.choose}}) when running interactively; in a non-interactive session (script, CI, \code{R CMD check}), \code{file} must be given explicitly.
#' @param dir The directory the file is downloaded (and, for a zip archive, extracted) into.
#' Defaults to the session's temporary directory, which is emptied when R exits. Pass an
#' explicit path (together with \code{cache = TRUE}) to keep the downloaded corpus between
#' sessions.
#' @param collapse Indicates whether or not lines of each documents should collapse together or not.
#' @param sep Separator between text fields.
#' @param categories Columns that should be considered as categorical data.
#' @param cache Whether the downloaded (and, for a zip archive, extracted) files are kept in
#' \code{dir} and reused on the next call. They are deleted, and downloaded again every time,
#' by default.
#' @return The text contained in the dowloaded file.
#' @export
#' @seealso \code{\link[utils]{download.file}}, \code{\link[utils]{unzip}}
#' @examples
#' # Not run automatically: this downloads a 31 MB archive from a third-party server, so it
#' # depends on both the network and that server staying up.
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' # Keep the archive between calls, in a directory of your choosing
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip", dir = "~/corpora", cache = TRUE)
#' }
loadtext <-
  function (file = NULL, dir = tempdir (), collapse = TRUE, sep = NULL, categories = NULL, cache = FALSE)
  {
    # Not the user's home directory: CRAN policy 1.6 reserves that for an explicit opt-in,
    # which 'dir' provides.
    dir = path.expand (dir)
    dir = sub ("(?<=.)/+$", "", dir, perl = TRUE)
    if (!dir.exists (dir))
      dir.create (dir, recursive = TRUE)
    if (is.null (file))
    {
      if (interactive ())
        file = file.choose ()
      else
        stop ("loadtext: 'file' must be specified in a non-interactive session (script, ",
              "R CMD check, CI, ...); the interactive file.choose() dialog is not available. ",
              "Please pass a file path or URL explicitly.")
    }
    mainfile = file
    download = grepl ("^https?://", file)
    if (download)
    {
      mainfile = file.path (dir, tail (strsplit (file, "/") [[1]], 1))
      if (!(cache && file.exists (mainfile)))
        utils::download.file (file, mainfile)
    }
    ext = tail (strsplit (mainfile, ".", fixed = TRUE) [[1]], 1)
    files = NULL
    if (ext %in% c ("zip"))
    {
      entries = utils::unzip (mainfile, exdir = dir, list = TRUE) [, 1]
      files = file.path (dir, entries)
      if (!(cache && all (file.exists (files))))
        utils::unzip (mainfile, exdir = dir, files = entries)
    }
    else
      files = mainfile
    corpus = NULL
    if (is.null (sep))
    {
      corpus = as.vector (sapply (files, function (file)
      {
        text = readLines (file, n = -1, warn = FALSE)
        if (collapse)
          text = paste (text, collapse = " ")
        return (text)
      }))
      corpus = corpus [!sapply (corpus, function (text) grepl ("^\\s*$", text))]
    }
    else
    {
      corpus = lapply (files, function (file)
      {
        text = utils::read.table (file, sep = sep, quote = "")
        return (text)
      })
      if (length (corpus) > 1)
        corpus = do.call (rbind, corpus)
      corpus [] = lapply(corpus, as.character)
      if (!is.null (categories))
        corpus [categories] = lapply(corpus [categories], factor)
    }

    if (download && !cache)
      file.remove (mainfile)
    if ((ext %in% c ("zip")) && !cache)
      sapply (files, function (file) file.remove (file))
    return (corpus)
  }

#' Plot word cloud
#'
#' Plot a word cloud based on the word frequencies in the documents.
#' @name plotcloud
#' @param corpus The corpus of documents (a vector of characters) or the vocabulary of the documents (result of function \code{getvocab}).
#' @param k A categorical variable (vector or factor).
#' @param stopwords The language whose stop words are removed (\code{"en"}, ...), or
#' \code{NULL} to keep them. A list of words of your own goes to \code{excludewords}.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{plotzipf}}, \code{\link{getvocab}}, \code{\link[wordcloud]{wordcloud}}
#' @examples
#' data (capitals)
#' plotcloud (capitals)
#' vocab = getvocab (capitals, mincount = 1, lang = NULL, stopwords = "en")
#' plotcloud (vocab)
plotcloud <-
  function (corpus, k = NULL, stopwords = "en", ...)
  {
    l = NULL
    labels = NULL
    kk = 1
    if (is.null (k))
      l = list (corpus)
    else
    {
      kk = sort (unique (k))
      for (i in kk)
        l = c (l, list (corpus [k == i]))
      if (is.factor (k))
        labels = levels (k)
      else
        labels = paste ("Cluster", kk)
    }
    n = length (kk)
    nrow = round (sqrt (n))
    ncol = ceiling (n / nrow)
    graphics::layout (matrix (1:(nrow * ncol), ncol = ncol, byrow = TRUE))
    on.exit (graphics::layout (1))
    for (i in 1:n)
    {
      vocab = NULL
      freq = NULL
      words = NULL
      if ("text2vec_vocabulary" %in% class (l [[i]]))
      {
        words = l [[i]] [, "term"]
        freq = l [[i]] [, "term_count"]
      }
      else
      {
        vocab = getvocab (l [[i]], mincount = 1, stopwords = stopwords, lang = NULL)
        words = vocab [, "term"]
        freq = vocab [, "term_count"]
      }
      maxfreq = max (freq)
      col = unique (grDevices::gray (1 - ((tail (freq, 200) + maxfreq) / (maxfreq * 2))))
      wordcloud::wordcloud (words = words, freq = freq, min.freq = 1, max.words = 200, random.order = FALSE, rot.per = 1 / 3, colors = col)
      graphics::title (main = labels [i])
    }
  }

#' Plot rank versus frequency
#'
#' Plot the frequency of words in a document agains the ranks of those words. It also plot the Zipf law.
#' @name plotzipf
#' @param corpus The corpus of documents (a vector of characters) or the vocabulary of the documents (result of function \code{getvocab}).
#' @export
#' @seealso \code{\link{plotcloud}}, \code{\link{getvocab}}
#' @examples
#' data (capitals)
#' plotzipf (capitals)
#' vocab = getvocab (capitals, mincount = 1, lang = NULL)
#' plotzipf (vocab)
plotzipf <-
  function (corpus)
  {
    freq = NULL
    if ("text2vec_vocabulary" %in% class (corpus))
      freq = corpus [, "term_count"]
    else
      freq = getvocab (corpus, mincount = 1, stopwords = NULL, lang = NULL) [, "term_count"]
    rank = 1:length (freq)
    freq = freq [rev (rank)]
    logd = data.frame (logrank = log2 (rank), logfreq = log2 (freq))
    model = stats::lm (logfreq ~ logrank, weights = freq, data = logd)
    # 'scipen' is raised so the log-log axes are labelled 10000 rather than 1e+04, but it is
    # a global option: leaving it set silently changed the way every subsequent print() in the
    # user's session formats numbers (CRAN policy 1.6 forbids it).
    old = options (scipen = freq [1])
    on.exit (options (old))
    graphics::plot (x = rank, y = freq, log = "xy", xlab = "Rank", ylab = "Frequency", t = "l")
    graphics::lines (rank, 2^model$coefficients [1] / rank^(-model$coefficients [2]), col = "red", lty = 2)
    graphics::legend ("topright", col = 1:2, legend = c ("Observations", "Zipf's law"), lty = 1:2, bty = "n")
  }

#' Model predictions
#'
#' This function predicts values based upon a model trained for text mining.
#' @name predict.textmining
#' @param object The classification model (of class \code{\link{textmining-class}}, created by \code{\link{TEXTMINING}}.
#' @param test The test set (a \code{data.frame})
#' @param fuzzy A boolean indicating whether fuzzy classification is used or not.
#' @return A vector of predicted values (\code{factor}).
#' @param ... Other parameters.
#' @export
#' @method predict textmining
#' @seealso \code{\link{TEXTMINING}}, \code{\link{textmining-class}}
#' @examples
#' \donttest{
#' require (text2vec)
#' # A small subset of movie_review is used here so this example runs quickly.
#' data ("movie_review")
#' d = movie_review [1:300, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' model = TEXTMINING (d$train.x, NB, labels = d$train.y, mincount = 10)
#' pred = predict (model, d$test.x)
#' evaluation (pred, d$test.y)
#' }
predict.textmining <- function (object, test, fuzzy = FALSE, ...)
{
  if (is.null (object$vectorizer))
    stop ("predict.textmining: this model was built with TEXTMINING (..., vector = \"words\"), ",
          "which vectorises the vocabulary rather than the documents. There is no document ",
          "vectorizer to project 'test' onto, so new documents cannot be predicted. Use ",
          "vector = \"docs\" if you need to classify unseen documents.")
  test = vectorize.docs (corpus = test, vectorizer = object$vectorizer)
  return (predict (object$res, as.matrix (test), fuzzy, ...))
}

#' Document query
#'
#' Search for documents similar to the query.
#' @name query.docs
#' @param docvectors The vectorized documents.
#' @param query The query (vectorized or raw text).
#' @param vectorizer The vectorizer that has been used to vectorize the documents.
#' @param nres The number of results.
#' @return The indices of the documents the most similar to the query.
#' @export
#' @seealso \code{\link{vectorize.docs}}, \code{\link[text2vec]{sim2}}
#' @examples
#' \donttest{
#' require (text2vec)
#' # A small subset of movie_review is used here so this example runs quickly.
#' data (movie_review)
#' reviews = movie_review$review [1:300]
#' vectorizer = vectorize.docs (corpus = reviews, returndata = FALSE)
#' docs = vectorize.docs (corpus = reviews, vectorizer = vectorizer)
#' query.docs (docs, reviews [1], vectorizer)
#' query.docs (docs, docs [1, ], vectorizer)
#' }
query.docs <-
  function (docvectors, query, vectorizer, nres = 5)
  {
    if (is.character (query))
      query = vectorize.docs (vectorizer, query)
    # vectorize.docs() (and a single-row slice of its result, e.g. docvectors [1, ]) returns a
    # data.frame; text2vec::sim2() requires proper matrix/Matrix objects for both 'x' and 'y',
    # so both need to be coerced explicitly (a bare matrix (query, nrow = 1) on a data.frame does
    # not flatten it to a numeric row, it produces a one-row list-matrix).
    docvectors = as.matrix (docvectors)
    query = matrix (as.numeric (as.matrix (query)), nrow = 1)
    taboo = apply (docvectors, 1, function (v) all (v == query))
    return (names (head (sort (text2vec::sim2 (x = docvectors [!taboo, ], y = query, method = "cosine", norm = "l2") [, 1], decreasing = TRUE), nres)))
  }

#' Word query
#'
#' Search for words similar to the query.
#' @name query.words
#' @param wordvectors The vectorized words
#' @param origin The query (character).
#' @param sub Words to be substrated to the origin.
#' @param add Words to be Added to the origin.
#' @param nres The number of results.
#' @param lang The language of the words (NULL if no stemming).
#' @return The Words the most similar to the query.
#' @export
#' @seealso \code{\link{vectorize.words}}, \code{\link[text2vec]{sim2}}
#' @examples
#' \donttest{
#' # 'capitals' is small, so the word vectors are coarse and 'ndim' is reduced
#' # accordingly; phrase detection needs a much larger corpus.
#' data (capitals)
#' words = vectorize.words (capitals, mincount = 2, ndim = 10, maxiter = 5)
#' query.words (words, origin = "paris", sub = "france", add = "germany")
#' query.words (words, origin = "berlin", sub = "germany", add = "france")
#' }
query.words <-
  function (wordvectors, origin, sub = NULL, add = NULL, nres = 5, lang = "en")
  {
    # vectorize.words() returns a data.frame; text2vec::sim2() requires proper matrix/Matrix
    # objects for both 'x' and 'y', so coerce once, up front, so that every row-slice and
    # arithmetic operation on wordvectors below (including the resulting query vector 'q')
    # stays a matrix rather than silently degrading back to a data.frame.
    wordvectors = as.matrix (wordvectors)
    words = rownames (wordvectors)
    origin = intersect (words, SnowballC::wordStem (tolower (origin), language = lang))
    if (length (origin) == 0)
      return (character (0))
    if (!is.null (sub))
      sub = intersect (words, SnowballC::wordStem (tolower (sub), language = lang))
    if (!is.null (add))
      add = intersect (words, SnowballC::wordStem (tolower (add), language = lang))
    taboo = which (words %in% c (origin, sub, add))
    q = wordvectors [origin [1], , drop = FALSE]
    if ((!is.null (sub)) && (length (sub) > 0))
      q = q - apply (wordvectors [sub, , drop = FALSE], 2, sum)
    if ((!is.null (add)) && (length (add) > 0))
      q = q + apply (wordvectors [add, , drop = FALSE], 2, sum)
    return (names (head (sort (text2vec::sim2 (x = wordvectors [-taboo, ], y = q, method = "cosine", norm = "l2") [, 1], decreasing = TRUE), nres)))
  }

#' @keywords internal
stemtokenizer <-
  function (x, lang = "en")
  {
    tokens = text2vec::word_tokenizer (x)
    res = lapply (tokens, SnowballC::wordStem, language = lang)
    return (res)
  }

#' Text mining
#'
#' Apply data mining function on vectorized text
#' @name TEXTMINING
#' @param corpus The corpus.
#' @param miningmethod The data mining method.
#' @param vector Indicates the type of vectorization, documents (TF-IDF) or words (GloVe).
#' @param ... Parameters passed to the vectorisation and to the data mining method.
#' @return The result of the data mining method.
#' @export
#' @seealso \code{\link{predict.textmining}}, \code{\link{textmining-class}}, \code{\link{vectorize.docs}}, \code{\link{vectorize.words}}
#' @examples
#' \donttest{
#' require (text2vec)
#' # A small subset of movie_review is used here so this example runs quickly.
#' data ("movie_review")
#' d = movie_review [1:300, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' model = TEXTMINING (d$train.x, NB, labels = d$train.y, mincount = 10)
#' pred = predict (model, d$test.x)
#' evaluation (pred, d$test.y)
#' data (capitals)
#' clusters = TEXTMINING (capitals, HCA, vector = "words", k = 5, mincount = 2, ndim = 10, maxiter = 5)
#' plotclus (clusters$res, capitals, type = "tree", labels = TRUE)
#' }
TEXTMINING <-
  function (corpus, miningmethod, vector = c ("docs", "words"), ...)
  {
    if (vector [1] == "docs")
    {
      vectorizer = vectorize.docs (corpus = corpus, returndata = FALSE, ...)
      d = as.matrix (vectorize.docs (corpus = corpus, vectorizer = vectorizer))
      res = miningmethod (d, ...)
      res = list (vectorizer = vectorizer, vectors = d, res = res, vector = "docs")
      class (res) = "textmining"
    }
    else
    {
      d = as.matrix (vectorize.words (corpus = corpus, ...))
      res = miningmethod (d, ...)
      # No 'vectorizer' slot here: the "words" mode vectorises the vocabulary, not the
      # documents, so there is nothing to project a new document onto. The mode is recorded
      # so that predict.textmining() can say so instead of failing on a NULL vectorizer.
      res = list (vectors = d, res = res, vector = "words")
      class (res) = "textmining"
    }
    return (res)
  }

#' @keywords internal
tokens <-
  function (corpus, lang = NULL, removesinglechars = TRUE)
  {
    tokenizer = text2vec::word_tokenizer
    if (!is.null (lang))
      tokenizer = stemtokenizer
    preprocessor = function (text) cleanup (text, removesinglechars = removesinglechars)
    return (text2vec::itoken (corpus, preprocessor = preprocessor, tokenizer = tokenizer, ids = 1:length (corpus), progressbar = FALSE, lang = lang))
  }

#' Document vectorization
#'
#' Vectorize a corpus of documents.
#' @name vectorize.docs
#' @param vectorizer The document vectorizer.
#' @inheritParams getvocab
#' @param transform Transformation (TF-IDF, LSA, L1 normanization, or nothing).
#' @param latentdim Number of latent dimensions if LSA transformation is performed.
#' @param returndata If true, the vectorized documents are returned. If false, a "vectorizer" is returned.
#' @param sparse Whether the document-term matrix is returned as a sparse matrix
#' (\code{dgCMatrix}) rather than as an ordinary \code{data.frame}. A document-term matrix is
#' mostly zeros, and storing them all takes about 4.5 GB for 20000 documents and 30000 terms,
#' so anything but a small corpus needs \code{TRUE}. Every method of the package accepts
#' either.
#' @param ... Other parameters.
#' @return The vectorized documents, as a \code{data.frame} or, if \code{sparse} is
#' \code{TRUE}, as a sparse matrix.
#' @export
#' @seealso \code{\link{query.docs}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{vectorizers}}
#' @examples
#' \donttest{
#' require (text2vec)
#' # A small subset of movie_review is used here so this example runs quickly.
#' data ("movie_review")
#' reviews = movie_review [1:300, ]
#' # Clustering
#' docs = vectorize.docs (corpus = reviews$review, transform = "tfidf")
#' km = KMEANS (docs [sample (nrow (docs), 50), ], k = 10)
#' # Classification
#' d = reviews [, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' vectorizer = vectorize.docs (corpus = d$train.x,
#'                              returndata = FALSE, mincount = 10)
#' train = vectorize.docs (corpus = d$train.x, vectorizer = vectorizer)
#' test = vectorize.docs (corpus = d$test.x, vectorizer = vectorizer)
#' model = NB (as.matrix (train), d$train.y)
#' pred = predict (model, as.matrix (test))
#' evaluation (pred, d$test.y)
#' }
vectorize.docs <-
  function (vectorizer = NULL, corpus = NULL, lang = "en", stopwords = lang, excludewords = NULL, ngram = 1, mincount = 10, minphrasecount = NULL, transform = c ("tfidf", "lsa", "l1", "none"), latentdim = 50, returndata = TRUE, removesinglechars = TRUE, sparse = FALSE, ...)
  {
    if (is.null (vectorizer))
      vectorizer = createvectorizer (corpus, lang = lang, stopwords = stopwords, excludewords = excludewords, ngram = ngram, mincount = mincount, minphrasecount = minphrasecount, transform = transform, latentdim = latentdim, removesinglechars = removesinglechars)
    if (returndata)
    {
      it = NULL
      if (is.null (corpus))
        it = vectorizer$tokens
      else
      {
        it = tokens (corpus, lang = lang, removesinglechars = removesinglechars)
        if (!is.null (vectorizer$phrases))
          it = vectorizer$phrases$transform (it)
      }
      dtm = text2vec::create_dtm (it, vectorizer$vectorizer)
      if (vectorizer$transform == "l1")
        dtm = text2vec::normalize (dtm, "l1")
      else if (vectorizer$transform == "tfidf")
        dtm = vectorizer$tfidf$transform(dtm)
      else if (vectorizer$transform == "lsa")
      {
        dtm = vectorizer$tfidf$transform(dtm)
        dtm = vectorizer$lsa$transform(dtm)
      }
      # A document-term matrix is naturally sparse, and text2vec produces it as such.
      # as.data.frame (as.matrix (.)) fills in every zero: on a realistic corpus of 20000
      # documents and 30000 terms that is a 4.5 GB data.frame. Keeping the default as it was
      # (students index the result like an ordinary data.frame), but sparse = TRUE returns
      # the sparse matrix, which every method of the package that takes a matrix accepts.
      if (sparse)
        return (dtm)
      return (as.data.frame (as.matrix (dtm)))
    }
    else
      return (vectorizer)
  }

#' Word vectorization
#'
#' Vectorize words from a corpus of documents.
#' @name vectorize.words
#' @param ndim The number of dimensions of the vector space.
#' @param maxwords The maximum number of words.
#' @inheritParams getvocab
#' @param window Window for term-co-occurrence matrix construction.
#' @param maxcooc Maximum number of co-occurrences to use in the weighting function.
#' @param maxiter The maximum number of iteration to fit the GloVe model.
#' @param epsilon Defines early stopping strategy when fit the GloVe model.
#' @param ... Other parameters.
#' @return The vectorized words.
#' @export
#' @seealso \code{\link{query.words}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{vectorizers}}
#' @examples
#' \donttest{
#' # 'capitals' is small, so the word vectors are coarse and 'ndim' is reduced
#' # accordingly; phrase detection needs a much larger corpus.
#' data (capitals)
#' words = vectorize.words (capitals, mincount = 2, ndim = 10, maxiter = 5)
#' query.words (words, origin = "paris", sub = "france", add = "germany")
#' query.words (words, origin = "berlin", sub = "germany", add = "france")
#' }
vectorize.words <-
  function (corpus = NULL, ndim = 50, maxwords = NULL, mincount = 5, minphrasecount = NULL, window = 5, maxcooc = 10, maxiter = 10, epsilon = 0.01, lang = "en", stopwords = lang, excludewords = NULL, removesinglechars = TRUE, ...)
  {
    it = createiterator (corpus, lang = lang, removesinglechars = removesinglechars)
    phrases = NULL
    if ((!is.null (minphrasecount)) && (minphrasecount > 0))
    {
      phrases = addphrases (it, mincount = minphrasecount)
      it = phrases$transform (it)
    }
    vocab = getvocab (corpus, mincount = mincount, minphrasecount = minphrasecount, ngram = 1, stopwords = stopwords, excludewords = excludewords, it = it, lang = lang)
    vectorizer = createvectorizer (corpus, it = it, phrases = phrases, vocab = vocab, stopwords = stopwords, excludewords = excludewords, ngram = 1, mincount = mincount, minphrasecount = minphrasecount, removesinglechars = removesinglechars)
    tcm = text2vec::create_tcm (vectorizer$tokens, vectorizer$vectorizer, skip_grams_window = window)
    glove = text2vec::GlobalVectors$new (rank = ndim, x_max = maxcooc)
    words = glove$fit_transform (tcm, n_iter = maxiter, convergence_tol = epsilon)
    words = words + t (glove$components)
    if (!is.null (maxwords))
    {
      fw = frequentwords (vocab, maxwords)
      words = words [fw, ]
    }
    return (as.data.frame (words))
  }
