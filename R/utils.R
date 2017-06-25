#' cat0, cat with sep="", similar to paste0
#' @param \ldots arguments passed to cat
#' @noRd
cat0 <- function(...) {
  cat(..., sep="", append=TRUE)
}

#' Create mass label
#'
#' Identifying the experiments by the running time/order is complicated.
#' Sometimes the instrument records a new run with the same settings which moves
#' all indicies. Same is true for the start times.
#'
#' We modify the target mass slightly by rounding to one decimal place and use
#' the second to fourth (default) to encode the id. Subsequently it is possible
#' to find the results by their individual mass.
#'
#' @param x double, original mass
#' @param id double, run id
#' @param divisor double, divisor (determines which decimal place)
#' @return double, mass label (id encoded in the second to fourth decimal place)
#' @seealso \code{\link{.massLabelToId}}
#' @noRd
#' @example
#' library("topdown")
#' topdown:::.massLabel(c(750, 1000), c(1, 100))
.massLabel <- function(x, id, divisor=10000L) {
  if (any(log10(divisor) <= log10(id) + 1L)) {
    stop(sQuote("divisor"), " has to be at least two digits more than ",
         sQuote("id"))
  }
  round(x, 1L) + id / divisor
}

#' Extract ID from mass labels
#'
#' @param x double, mass label
#' @param multiplicator double, (determines which decimal place)
#' @seealso \code{\link{.massLabel}}
#' @noRd
.massLabelToId <- function(x, multiplicator=10000L) {
  x <- as.double(x)
  (x - round(x, 1)) * multiplicator
}

#' verbose output
#' @param \ldots arguments passed to message
#' @noRd
.msg <- function(verbose, ...) {
  if (verbose) {
    message(...)
  }
}

#' similar to lengths but for rows
#' @param x list of data.frames/matrices
#' @noRd
.nrows <- function(x) {
  stopifnot(is.list(x))
  vapply(x, nrow, double(1L))
}

#' The ScanHeadsMan output for the scan conditons contains a column
#' TargetedMassList with the format "(mz=[0-9]+\.[0-9]+ z=[0-9]+ name=)". This
#' function converts this format to truncated (one decimal place) mz values.
#' @param x character
#' @return double
#' @noRd
.targetedMassListToMz <- function(x) {
  stopifnot(is.character(x))
  trunc(as.double(gsub("^.*mz=([^ ]+) z.*$", "\\1", x)) * 10L) / 10L
}
