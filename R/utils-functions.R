#' cat0, cat with sep="", similar to paste0
#' @param \ldots arguments passed to cat
#' @noRd
cat0 <- function(...) {
  cat(..., sep="", append=TRUE)
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
