#' @keywords internal
addalpha <-
  function (colors, a = 32)
  {
    return (sapply (colors, function (color)
    {
      rgb = grDevices::col2rgb (color)
      r = rgb [1]
      g = rgb [2]
      b = rgb [3]
      result = c (result, grDevices::rgb (r, g, b, a, maxColorValue = 255))
    }))
  }

#' Running time
#'
#' Return the running time of a function
#' @name runningtime
#' @param FUN The function to be evaluated.
#' @param ... The parameters to be passes to function \code{FUN}.
#' @return The running time of function \code{FUN}.
#' @export
#' @seealso \code{\link[base]{difftime}}
#' @examples
#' sqrt (x = 1:100)
#' runningtime (sqrt, x = 1:100)
runningtime <-
  function (FUN, ...)
  {
    start = Sys.time ()
    FUN (...)
    end = Sys.time ()
    return (end - start)
  }

#' Splits a dataset into training set and test set
#'
#' This function splits a dataset into training set and test set. Return an object of class \code{\link{dataset-class}}.
#' @name splitdata
#' @param dataset The dataset to be split (\code{data.frame} or \code{matrix}).
#' @param target The column index of the target variable (class label or response variable).
#' @param size The size of the training set (as an integer value).
#' @param seed A specified seed for random number generation.
#' @return An object of class \code{\link{dataset-class}}.
#' @export
#' @seealso \code{\link{dataset-class}}
#' @examples
#' require (datasets)
#' data (iris)
#' d = splitdata (iris, 5)
#' str (d)
splitdata <-
  function (dataset, target, size = round (0.7 * nrow (dataset)), seed = NULL)
  {
    set.seed (seed)
    s = sample (nrow (dataset), size)
    train.x = dataset [s, -target]
    train.y = dataset [s, target]
    test.x = dataset [-s, -target]
    test.y = dataset [-s, target]
    res = list (train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y)
    class (res) = "dataset"
    return (res)
  }
