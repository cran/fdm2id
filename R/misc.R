fdm2id.globals <- new.env (emptyenv ())
fdm2id.globals$export <- TRUE

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

#' Close a graphics device
#'
#' Close the graphics device driver
#' @name closegraphics
#' @seealso \code{\link{exportgraphics}}, \code{\link{toggleexport}}, \code{\link[grDevices]{dev.off}}
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' }
closegraphics <-
  function ()
  {
    if (fdm2id.globals$export)
      grDevices::dev.off ()
  }

#' Open a graphics device
#'
#' Starts the graphics device driver
#' @name exportgraphics
#' @param file A character string giving the name of the file.
#' @param type The type of graphics device.
#' @param ... Other parameters.
#' @seealso \code{\link{closegraphics}}, \code{\link{toggleexport}}, \code{\link[grDevices]{Devices}}
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' }
exportgraphics <-
  function (file, type = tail (strsplit (file, split = "\\.") [[1]], 1), ...)
  {
    if (is.character (type))
      type = get (type)
    if (fdm2id.globals$export)
      type (file, ...)
  }

#' @rdname toggleexport
#' @export
exportgraphics.off <-
  function ()
  {
    toggleexport (FALSE)
  }

#' @rdname toggleexport
#' @export
exportgraphics.on <-
  function ()
  {
    toggleexport (TRUE)
  }

#' Rotation
#'
#' Rotation on two variables of a numeric dataset
#' @name rotation
#' @param d The dataset.
#' @param angle The angle of the rotation.
#' @param axis The axis.
#' @param range The range of the angle (360, 2*pi, 100, ...)
#' @export
#' @examples
#' d = data.parabol ()
#' d [, -3] = rotation (d [, -3], 45, range = 360)
#' plotdata (d [, -3], d [, 3])
rotation = function (d, angle, axis = 1:2, range = 2 * pi)
{
  theta = 2 * pi * angle / range
  rot = diag (ncol (d))
  rot [axis, axis] = matrix (c (cos (theta), sin (theta), -sin (theta), cos (theta)), ncol = 2)
  res = as.matrix (d) %*% rot
  return (res)
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
    if (size < 1)
      size = round (size * nrow (dataset))
    s = sample (nrow (dataset), size)
    train.x = dataset [s, -target]
    train.y = dataset [s, target]
    test.x = dataset [-s, -target]
    test.y = dataset [-s, target]
    res = list (train.x = train.x, train.y = train.y, test.x = test.x, test.y = test.y)
    class (res) = "dataset"
    return (res)
  }

#' Toggle graphic exports
#'
#' Toggle graphic exports on and off
#' @name toggleexport
#' @aliases exportgraphics.off exportgraphics.on toggleexport.off toggleexport.on
#' @param export If \code{TRUE}, exports are activated, if \code{FALSE}, exports are deactivated. If \code{null}, switches on and off.
#' @seealso \code{\link{closegraphics}}, \code{\link{exportgraphics}}
#' @rdname toggleexport
#' @export
#' @examples
#' \dontrun{
#' data (iris)
#' toggleexport (FALSE)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' toggleexport (TRUE)
#' exportgraphics ("export.pdf")
#' plotdata (iris [, -5], iris [, 5])
#' closegraphics()
#' }
toggleexport <-
  function (export = NULL)
  {
    if (is.null (export))
      export = !fdm2id.globals$export
    fdm2id.globals$export = export
  }

#' @rdname toggleexport
#' @export
toggleexport.off <-
  function ()
  {
    toggleexport (FALSE)
  }

#' @rdname toggleexport
#' @export
toggleexport.on <-
  function ()
  {
    toggleexport (TRUE)
  }
