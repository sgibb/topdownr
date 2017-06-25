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

#' The ScanHeadsMan output for the scan conditons contains a column
#' TargetedMassList with the format "(mz=[0-9]+\.[0-9]+ z=[0-9]+ name=)". This
#' function converts this format to truncated (one decimal place) mz values.
#' @param x character
#' @return double
#' @noRd
.targetMassList2Mz <- function(x) {
  stopifnot(is.character(x))
  trunc(as.double(gsub("^.*mz=([^ ]+) z.*$", "\\1", x)) * 10L) / 10L
}
