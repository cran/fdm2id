#' Document vectorization object
#'
#' This class contains a vectorization model for textual documents.
#' @name vectorizer-class
#' @slot vectorizer The vectorizer.
#' @slot transform The transformation to be applied after vectorization (normalization, TF-IDF).
#' @slot phrases The phrase detection method.
#' @slot tfidf The TF-IDF transformation.
#' @slot lsa The LSA transformation.
#' @slot tokens The token from the original document.
#' @exportClass vectorizer
#' @seealso \code{\link{vectorize.docs}}, \code{\link{query.docs}}
setClass ("vectorizer",
          representation (vectorizer = "function",
                          transform = "character",
                          phrases = "ANY",
                          tfidf = "ANY",
                          lsa = "ANY",
                          tokens = "ANY"))

#' Text mining object
#'
#' Object used for text mining.
#' @name textmining-class
#' @slot vectorizer The vectorizer.
#' @slot vectors The vectorized dataset.
#' @slot res The result of the text mining method.
#' @exportClass textmining
#' @seealso \code{\link{TEXTMINING}}, \code{\link{vectorize.docs}}
setClass ("textmining",
          representation (vectorizer = "function",
                          vectors = "matrix",
                          res = "ANY"))

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
  function (corpus)
  {
    res = sapply (corpus, function (text) tolower (text))
    res = sapply (res, function (text) gsub ("[^[:alnum:]]", " ", text))
    res = sapply (res, function (text) gsub ("\\b[[:alnum:]]{1}\\b", "", text))
    res = sapply (res, function (text) gsub ("\\s+", " ", text))
    return (res)
  }

#' @keywords internal
createiterator <-
  function (corpus, lang,  minphrasecount = NULL)
  {
    it = tokens (corpus, lang = lang)
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
  function (corpus, it = NULL, phrases = NULL, vocab = NULL, lang, stopwords = lang, ngram = 1, mincount = 10, minphrasecount = NULL,
            transform = c ("none", "l1", "tfidf", "lsa"), latentdim = 50)
  {
    if (is.null (it))
      it = createiterator (corpus, lang, minphrasecount)
    if (is.null (vocab))
      vocab = getvocab (corpus, mincount, minphrasecount, ngram, stopwords, it = it, lang = lang)
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
#' @param mincount Minimum word count to be considered as frequent.
#' @param minphrasecount Minimum collocation of words count to be considered as frequent.
#' @param ngram maximum size of n-grams.
#' @param lang The language of the documents (NULL if no stemming).
#' @param stopwords Stopwords, or the language of the documents. NULL if stop words should not be removed.
#' @return The most frequent words of the corpus.
#' @export
#' @seealso \code{\link{getvocab}}
#' @examples
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' frequentwords (text, 100)
#' vocab = getvocab (text)
#' frequentwords (vocab, 100)
#' }
frequentwords <-
  function (corpus, nb, mincount = 5, minphrasecount = NULL, ngram = 1, lang = "en", stopwords = lang)
  {
    vocab = NULL
    if ("text2vec_vocabulary" %in% class (corpus))
      vocab = corpus
    else
      vocab = getvocab (corpus, mincount = mincount, minphrasecount = minphrasecount, ngram = ngram, lang = lang, stopwords = stopwords)
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
#' @param stopwords Stopwords, or the language of the documents. NULL if stop words should not be removed.
#' @param ... Other parameters.
#' @return The vocabulary used in the corpus of documents.
#' @export
#' @seealso \code{\link{plotzipf}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{create_vocabulary}}
#' @examples
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' vocab1 = getvocab (text) # With stemming
#' nrow (vocab1)
#' vocab2 = getvocab (text, lang = NULL) # Without stemming
#' nrow (vocab2)
#' }
getvocab <-
  function (corpus, mincount = 5, minphrasecount = NULL, ngram = 1, lang = "en", stopwords = lang, ...)
  {
    dots = list (...)
    it = NULL
    if (!is.null (dots$it))
      it = dots$it
    else
    {
      it = tokens (corpus, lang = lang)
      if ((!is.null (minphrasecount)) && (minphrasecount > 0))
      {
        phrases = addphrases (it, mincount = minphrasecount)
        it = phrases$transform (it)
      }
    }
    sw = character(0)
    if (!is.null (stopwords))
    {
      if (length (stopwords) == 1)
        sw = stopwords::stopwords (stopwords)
      else
        sw = stopwords
    }
    vocab = text2vec::create_vocabulary (it, ngram = c (1, ngram), stopwords = sw)
    vocab = text2vec::prune_vocabulary (vocab, term_count_min = mincount)
    return (vocab)
  }

#' load a text file
#'
#' (Down)Load a text file (and extract it if it is in a zip file).
#' @name loadtext
#' @param file The path or URL of the text file.
#' @param dir The (temporary) directory, where the file is downloaded. The file is deleted at the end of this function.
#' @param collapse Indicates whether or not lines of each documents should collapse together or not.
#' @return The text contained in the dowloaded file.
#' @export
#' @seealso \code{\link[utils]{download.file}}, \code{\link[utils]{unzip}}
#' @examples
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' }
loadtext <-
  function (file = file.choose (), dir = "~/", collapse = TRUE)
  {
    mainfile = file
    download = grepl ("^https?://", file)
    if (download)
    {
      mainfile = paste (dir, tail (strsplit (file, "/") [[1]], 1), sep = "")
      utils::download.file (file, mainfile)
    }
    ext = tail (strsplit (mainfile, ".", fixed = TRUE) [[1]], 1)
    files = NULL
    if (ext %in% c ("zip"))
    {
      files = utils::unzip (mainfile, exdir = dir, list = TRUE) [, 1]
      utils::unzip (mainfile, exdir = dir, files = files)
      files = paste (dir, files, sep = "")
    }
    else
      files = mainfile
    corpus = as.vector (sapply (files, function (file)
    {
      text = readLines (file, n = -1, warn = FALSE)
      if (collapse)
        text = paste (text, collapse = " ")
      return (text)
    }))
    corpus = corpus [!sapply (corpus, function (text) grepl ("^\\s*$", text))]
    if (download)
      file.remove (mainfile)
    if (ext %in% c ("zip"))
      sapply (files, function (file) file.remove (file))
    return (corpus)
  }

#' Plot word cloud
#'
#' Plot a word cloud based on the word frequencies in the documents.
#' @name plotcloud
#' @param corpus The corpus of documents (a vector of characters) or the vocabulary of the documents (result of function \code{getvocab}).
#' @param k A categorial variable (vector or factor).
#' @param stopwords Stopwords, or the language of the documents. NULL if stop words should not be removed.
#' @param ... Other parameters.
#' @export
#' @seealso \code{\link{plotzipf}}, \code{\link{getvocab}}, \code{\link[wordcloud]{wordcloud}}
#' @examples
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' plotcloud (text)
#' vocab = getvocab (text, mincount = 1, lang = NULL, stopwords = "en")
#' plotcloud (vocab)
#' }
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
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' plotzipf (text)
#' vocab = getvocab (text, mincount = 1, lang = NULL)
#' plotzipf (vocab)
#' }
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
    options (scipen = freq [1])
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
#' \dontrun{
#' require (text2vec)
#' data ("movie_review")
#' d = movie_review [, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' model = TEXTMINING (d$train.x, NB, labels = d$train.y, mincount = 50)
#' pred = predict (model, d$test.x)
#' evaluation (pred, d$test.y)
#' }
predict.textmining <- function (object, test, fuzzy = FALSE, ...)
{
  test = vectorize.docs (corpus = test, vectorizer = object$vectorizer)
  return (predict (object$res, as.matrix (test), fuzzy, ...))
}

#' Document query
#'
#' Search for documents similar to the query.
#' @name query.docs
#' @param docvectors The vectorized documents.
#' @param query The query (vectorized or raw text).
#' @param vectorizer The vectorizer taht has been used to vectorize the documents.
#' @param nres The number of results.
#' @return The indices of the documents the most similar to the query.
#' @export
#' @seealso \code{\link{vectorize.docs}}, \code{\link[text2vec]{sim2}}
#' @examples
#' \dontrun{
#' require (text2vec)
#' data (movie_review)
#' vectorizer = vectorize.docs (corpus = movie_review$review,
#'                              minphrasecount = 50, returndata = FALSE)
#' docs = vectorize.docs (corpus = movie_review$review, vectorizer = vectorizer)
#' query.docs (docs, movie_review$review [1], vectorizer)
#' query.docs (docs, docs [1, ], vectorizer)
#' }
query.docs <-
  function (docvectors, query, vectorizer, nres = 5)
  {
    if (is.character (query))
      query = vectorize.docs (vectorizer, query)
    taboo = apply (docvectors, 1, function (v) all (v == query))
    return (names (head (sort (text2vec::sim2 (x = docvectors [!taboo, ], y = matrix (query, nrow = 1), method = "cosine", norm = "l2") [, 1], decreasing = TRUE), nres)))
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
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' words = vectorize.words (text, minphrasecount = 50)
#' query.words (words, origin = "paris", sub = "france", add = "germany")
#' query.words (words, origin = "berlin", sub = "germany", add = "france")
#' query.words (words, origin = "new_zealand")
#' }
query.words <-
  function (wordvectors, origin, sub = NULL, add = NULL, nres = 5, lang = "en")
  {
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
#' \dontrun{
#' require (text2vec)
#' data ("movie_review")
#' d = movie_review [, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' model = TEXTMINING (d$train.x, NB, labels = d$train.y, mincount = 50)
#' pred = predict (model, d$test.x)
#' evaluation (pred, d$test.y)
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' clusters = TEXTMINING (text, HCA, vector = "words", k = 9, maxwords = 100)
#' plotclus (clusters$res, text, type = "tree", labels = TRUE)
#' }
TEXTMINING <-
  function (corpus, miningmethod, vector = c ("docs", "words"), ...)
  {
    if (vector [1] == "docs")
    {
      vectorizer = vectorize.docs (corpus = corpus, returndata = FALSE, ...)
      d = as.matrix (vectorize.docs (corpus = corpus, vectorizer = vectorizer))
      res = miningmethod (d, ...)
      res = list (vectorizer = vectorizer, vectors = d, res = res)
      class (res) = "textmining"
    }
    else
    {
      d = as.matrix (vectorize.words (corpus = corpus, ...))
      res = miningmethod (d, ...)
      res = list (vectors = d, res = res)
      class (res) = "textmining"
    }
    return (res)
  }

#' @keywords internal
tokens <-
  function (corpus, lang = NULL)
  {
    ids = NULL
    if (length (corpus) > 1)
      ids = 1:length (corpus)
    tokenizer = text2vec::word_tokenizer
    if (!is.null (lang))
      tokenizer = stemtokenizer
    return (text2vec::itoken (corpus, preprocessor = cleanup, tokenizer = tokenizer, ids = 1:length (corpus), progressbar = FALSE, lang = lang))
  }

#' Document vectorization
#'
#' Vectorize a corpus of documents.
#' @name vectorize.docs
#' @param vectorizer The document vectorizer.
#' @param corpus The corpus of documents (a vector of characters).
#' @param lang The language of the documents (NULL if no stemming).
#' @param stopwords Stopwords, or the language of the documents. NULL if stop words should not be removed.
#' @param ngram maximum size of n-grams.
#' @param mincount Minimum word count to be considered as frequent.
#' @param minphrasecount Minimum collocation of words count to be considered as frequent.
#' @param transform Transformation (TF-IDF, LSA, L1 normanization, or nothing).
#' @param latentdim Number of latent dimensions if LSA transformation is performed.
#' @param returndata If true, the vectorized documents are returned. If false, a "vectorizer" is returned.
#' @param ... Other parameters.
#' @return The vectorized documents.
#' @export
#' @seealso \code{\link{query.docs}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{vectorizers}}
#' @examples
#' \dontrun{
#' require (text2vec)
#' data ("movie_review")
#' # Clustering
#' docs = vectorize.docs (corpus = movie_review$review, transform = "tfidf")
#' km = KMEANS (docs [sample (nrow (docs), 100), ], k = 10)
#' # Classification
#' d = movie_review [, 2:3]
#' d [, 1] = factor (d [, 1])
#' d = splitdata (d, 1)
#' vectorizer = vectorize.docs (corpus = d$train.x,
#'                              returndata = FALSE, mincount = 50)
#' train = vectorize.docs (corpus = d$train.x, vectorizer = vectorizer)
#' test = vectorize.docs (corpus = d$test.x, vectorizer = vectorizer)
#' model = NB (as.matrix (train), d$train.y)
#' pred = predict (model, as.matrix (test))
#' evaluation (pred, d$test.y)
#' }
vectorize.docs <-
  function (vectorizer = NULL, corpus = NULL, lang = "en", stopwords = lang, ngram = 1, mincount = 10, minphrasecount = NULL, transform = c ("tfidf", "lsa", "l1", "none"), latentdim = 50, returndata = TRUE, ...)
  {
    if (is.null (vectorizer))
      vectorizer = createvectorizer (corpus, lang = lang, stopwords = stopwords, ngram = ngram, mincount = mincount, minphrasecount = minphrasecount, transform = transform, latentdim = latentdim)
    if (returndata)
    {
      it = NULL
      if (is.null (corpus))
        it = vectorizer$tokens
      else
      {
        it = tokens (corpus, lang = lang)
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
      return (dtm)
    }
    else
      return (vectorizer)
  }

#' Word vectorization
#'
#' Vectorize words from a corpus of documents.
#' @name vectorize.words
#' @param corpus The corpus of documents (a vector of characters).
#' @param ndim The number of dimensions of the vector space.
#' @param maxwords The maximum number of words.
#' @param mincount Minimum word count to be considered as frequent.
#' @param minphrasecount Minimum collocation of words count to be considered as frequent.
#' @param window Window for term-co-occurence matrix construction.
#' @param maxcooc Maximum number of co-occurrences to use in the weighting function.
#' @param maxiter The maximum number of iteration to fit the GloVe model.
#' @param epsilon Defines early stopping strategy when fit the GloVe model.
#' @param lang The language of the documents (NULL if no stemming).
#' @param stopwords Stopwords, or the language of the documents. NULL if stop words should not be removed.
#' @param ... Other parameters.
#' @return The vectorized words.
#' @export
#' @seealso \code{\link{query.words}}, \code{\link[stopwords]{stopwords}}, \code{\link[text2vec]{vectorizers}}
#' @examples
#' \dontrun{
#' text = loadtext ("http://mattmahoney.net/dc/text8.zip")
#' words = vectorize.words (text, minphrasecount = 50)
#' query.words (words, origin = "paris", sub = "france", add = "germany")
#' query.words (words, origin = "berlin", sub = "germany", add = "france")
#' query.words (words, origin = "new_zealand")
#' }
vectorize.words <-
  function (corpus = NULL, ndim = 50, maxwords = NULL, mincount = 5, minphrasecount = NULL, window = 5, maxcooc = 10, maxiter = 10, epsilon = 0.01, lang = "en", stopwords = lang, ...)
  {
    it = createiterator (corpus, lang = lang)
    phrases = NULL
    if ((!is.null (minphrasecount)) && (minphrasecount > 0))
    {
      phrases = addphrases (it, mincount = minphrasecount)
      it = phrases$transform (it)
    }
    vocab = getvocab (corpus, mincount = mincount, minphrasecount = minphrasecount, ngram = 1, stopwords = stopwords, it = it, lang = lang)
    vectorizer = createvectorizer (corpus, it = it, phrases = phrases, vocab = vocab, stopwords = stopwords, ngram = 1, mincount = mincount, minphrasecount = minphrasecount)
    tcm = text2vec::create_tcm (vectorizer$tokens, vectorizer$vectorizer, skip_grams_window = window)
    glove = text2vec::GlobalVectors$new (word_vectors_size = ndim, vocabulary = vocab, x_max = maxcooc)
    words = glove$fit_transform (tcm, n_iter = maxiter, convergence_tol = epsilon)
    words = words + t (glove$components)
    if (!is.null (maxwords))
    {
      fw = frequentwords (vocab, maxwords)
      words = words [fw, ]
    }
    return (words)
  }
