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
  .massLabelToId(gsub("^.*ms2 ([^@]+)\\@.*$", "\\1", x))
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

#' Get fragmentation method from {ETD,CID,HCD}Activation
#' @param x data.frame/matrix, 3 columns (ETD,CID,HCD)
#' @return vector, fragmentation method
#' @noRd
.fragmentationMethod <- function(x) {
  methods <- c("None", "ETD", "CID", "ETcid", "HCD", "EThcd", "HCD/CID", "All")
  v <- c(ETDActivation=1L, CIDActivation=2L, HCDActivation=4L)
  stopifnot(all(colnames(x) %in% names(v)))
  x <- x[, names(v)]
  apply(x, 1L, function(i)methods[sum(v[as.logical(i)]) + 1L])
}

#' Split list/data.frame
#' @param x data.frame, e.g. from .ms2Experiments
#' @param cols character, colnames used to split
#' @return list
#' @noRd
.groupBy <- function(x, cols) {
  split(x, .groupByLabels(x, cols))
}

#' Create group labels from data.frame columns
#' @param x data.frame, e.g. from .ms2Experiments
#' @param cols character, colnames used to split
#' @return list
#' @noRd
.groupByLabels <- function(x, cols) {
  if (length(cols) > 1L) {
    ## `interaction` doesn't handle NA values, so use `paste` instead
    do.call(paste, c(x[, cols], sep=":"))
  } else {
    as.character(x[, cols])
  }
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
#' @param x character, mass label
#' @param idDigits integer, number of digits behind the decimal place that
#' store id information (not mass information) from the last one (e.g. 3 if
#' the id is 123 and the mz is 900.0123)
#' @seealso \code{\link{.massLabel}}
#' @noRd
.massLabelToId <- function(x, idDigits=3L) {
  # was the following before, but results in round errors ("7" becomes 6L)
  # x <- as.double(x)
  # as.integer((x - round(x, 1L)) * multiplicator)
  n <- nchar(x)
  as.integer(substring(x, n - idDigits + 1L, n))
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

#' shortend string to width and place "..." in the middle
#' @param x character
#' @param width number of letters allowed/width of terminal
#' @return character
#' @noRd
.snippet <- function(x, width=getOption("width")) {
  nc <- nchar(x)
  w <- (width - 2L:3L) %/% 2L
  ifelse(nc <= width, x, paste0(substring(x, 1L, w[1L]), "...",
                                substring(x, nc - w[2L] + 1L, nc)))
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
#' @return character, regexp for file extensions
#' @noRd
.topDownFileExtRx <- function(type=c("cfmt", "csv", "fasta", "mzml", "txt",
                                     "raw", "all")) {
  type <- match.arg(type)
  ext <- c(csv="experiments\\.csv", fasta="fasta", mzml="mz[Mm][Ll]",
           raw="raw", txt="txt")
  sel <- switch(type,
                "all" = seq_along(ext),
                "cfmt" = c("csv", "fasta", "mzml", "txt"),
                type)
  paste0("\\.", ext[sel], "$", collapse="|")
}

#' Update assignmentTable MzId values.
#' Because of the initial removal of non-matching peaks, each subsetting will
#' result in an mz index 1:n. This function regenerates the MzId column.
#' @param x assignmentTable
#' @return data.table, updated MzId column
#' @noRd
.updateAssignmentTableMzId <- function(x) {
  x[, MzId:=seq_len(.N), by=SpectrumId]
  setkey(x, SpectrumId, FragmentId, MzId)
  x
}

#' wrapper around vapply for FUN.VALUE=double(1L)
#' @noRd
.vapply1d <- function(X, FUN, ..., USE.NAMES=FALSE) {
  vapply(X=X, FUN=FUN, FUN.VALUE=NA_real_, ..., USE.NAMES=USE.NAMES)
}

#' wrapper around vapply for FUN.VALUE=logical(1L)
#' @noRd
.vapply1l <- function(X, FUN, ..., USE.NAMES=FALSE) {
  vapply(X=X, FUN=FUN, FUN.VALUE=NA, ..., USE.NAMES=USE.NAMES)
}
