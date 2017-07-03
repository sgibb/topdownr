#' cat0, cat with sep="", similar to paste0
#' @param \ldots arguments passed to cat
#' @noRd
cat0 <- function(...) {
  cat(..., sep="", append=TRUE)
}

#' The ScanHeadsMan output for the header information contains a column
#' FilterString with the format "FTMS + p NSI Full ms2 [0-9]+\.[0-9]+@hcd35.00
#' [xxx-yyy]". This function converts this format to the ID stored in the mass
#' label.
#' @param x character
#' @return double
#' @noRd
.filterStringToId <- function(x) {
  stopifnot(is.character(x))
  .massLabelToId(as.double(gsub("^.*ms2 ([^@]+)\\@.*$", "\\1", x)))
}

#' Create (nearly) CamelCase names. Could not correct "AGC" to "Agc".
#' @param x character
#' @return character
#' @noRd
.formatNames <- function(x) {
  x <- gsub("[[:punct:]]+", "", x)
  unlist(lapply(strsplit(x, " "), function(s) {
    if (length(s) == 1L) {
      # don't convert already camelcased strings
      paste0(toupper(substring(s, 1L, 1L)), substring(s, 2L), collapse="")
    } else {
      paste0(toupper(substring(s, 1L, 1L)), tolower(substring(s, 2L)),
             collapse="")
    }
  }))
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
  .vapply1d(x, nrow)
}

#' swap file extensions
#' @param x character, file name
#' @param ext character, new extension
#' @return character
#' @noRd
.swapFileExt <- function(x, ext="meth") {
  paste(file_path_sans_ext(x), ext, sep=".")
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

#' TopDown file extensions.
#' @param type character, which file ext
#' @return regexp for file extensions
#' @noRd
.topDownFileExtRx <- function(type=c("cmt", "csv", "mzml", "txt", "raw",
                                     "all")) {
  type <- match.arg(type)
  ext <- c(csv="experiments\\.csv", mzml="mz[Mm][Ll]", raw="raw", txt="txt")
  sel <- switch(type,
                "all" = seq_along(ext),
                "cmt" = c("csv", "mzml", "txt"),
                type)
  paste0("\\.", ext[sel], "$", collapse="|")
}

#' wrapper around vapply for FUN.VALUE=double(1L)
#' @noRd
.vapply1d <- function(X, FUN, ..., USE.NAMES=FALSE) {
  vapply(X=X, FUN=FUN, FUN.VALUE=NA_real_, ..., USE.NAMES=USE.NAMES)
}
