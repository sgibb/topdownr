#' @param x TopDownExperiment
#' @noRd
setMethod("[", "TopDownExperiment", function(x, i, j="missing", ...,
                                             drop="missing") {

  if (missing(i)) {
    i <- seq_along(x)
  }
  if (!(is.logical(i) | is.numeric(i) | is.character(i))) {
    stop("subsetting works only with numeric, logical or character")
  }
  if (length(x) < length(i)) {
    stop("subscript out of bounds")
  }
  if (is.character(i)) {
    i <- unique(match(i, featureNames(x)))
  }
  if (anyNA(i)) {
    stop("subscript out of bounds")
  }
  if (is.numeric(i)) {
    if (max(i) > length(x) | min(i) < 1) {
      stop("subscript out of bounds")
    }
  }
  d0 <- dim(x)
  fn <- featureNames(x)[i]
  x@assignmentTable <- assignmentTable(x)[SpectrumId %in% fn,]

  ## overwrite MSnbase log msg (and duplicated log by TDE subsetting if j is
  ## given)
  processing <- x@processingData@processing
  x <- callNextMethod(x=x, i=i, j=j, ..., drop=drop)

  if (!missing(j) && length(j) && j[1L] != "missing") {
    if (is.character(j)) {
      x <- .filterFragmentIonOrType(x, j)
    } else if (is.numeric(j)) {
      x <- .filterFragmentPos(x, j)
    } else if (is.logical(j) && length(j) == nrow(fragmentTable(x))) {
      x <- .filterFragmentId(x, .fragmentId(x)[j])
    } else {
      stop("subscript out of bounds")
    }
  }

  d1 <- dim(x)
  x@processingData@processing <- processing
  x <- .logmsg(x, paste0("Subset [", d0[1L], ";", d0[2L], "] to [",
                                     d1[1L], ";", d1[2L], "]."))

  x
})

#' Accessor for assignmentTable, not exported yet
#' @param object TopDownExperiment
#' @return data.table
#' @noRd
setMethod("assignmentTable", "TopDownExperiment", function(object) {
  object@assignmentTable
})

#' Dim (number of spectra, number of spectra-fragment-assignments)
#' @param object TopDownExperiment
#' @return numeric
#' @noRd
setMethod("dim", "TopDownExperiment", function(x) {
  as.integer(c(length(x), nrow(assignmentTable(x))))
})

#' Accessor for fragmentTable, not exported yet
#' @param object TopDownExperiment
#' @return data.table
#' @noRd
setMethod("fragmentTable", "TopDownExperiment", function(object) {
  object@fragmentTable
})

#' @param object TopDownExperiment
#' @noRd
setMethod("show", "TopDownExperiment", function(object) {
  cat("TopDown Experiment data (\"", class(object), "\").\n", sep="")

  if (nchar(object@sequence)) {
    cat("- - - Protein data - - -\n")
    prefix <- sprintf("Amino acid sequence (%d):", nchar(object@sequence))
    cat(prefix, .snippet(object@sequence, getOption("width") - nchar(prefix)),
        "\n")
  }

  if (nrow(fragmentTable(object))) {
    cat("- - - Fragment data - - -\n")
    cat("Number of theoretical fragments:", nrow(fragmentTable(object)), "\n")
    fragments <- sort(unique(.fragmentTypes(object)))
    prefix <- sprintf("Theoretical fragment types (%d):", length(fragments))
    cat(prefix, .snippet(paste0(fragments, collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (nrow(assignmentTable(object))) {
    cat("- - - Spectra data - - -\n")
    cat("Number of Spectra:", length(object), "\n")
    cat("Number of Spectra-Fragment-Assignments:",
        nrow(assignmentTable(object)), "\n")
    fragments <- sort(unique(.fragmentTypes(object)[
      .fragmentId(object) %in% assignmentTable(object)$FragmentId]))
    prefix <- sprintf("Assigned fragment types (%d):", length(fragments))
    cat(prefix, .snippet(paste0(fragments, collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (length(processingData(object)@processing)) {
    cat("- - - Processing information - - -\n")
    cat(paste0(processingData(object)@processing, collapse="\n"), "\n")
  }

  invisible(NULL)
})

