#' @param x TopDownExperiment
#' @noRd
setMethod("[", "TopDownExperiment", function(x, i, j="missing", ...,
                                             drop="missing") {

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
  fn <- featureNames(x)[i]
  x@assignmentTable <- assignmentTable(x)[SpectrumId %in% fn,]

  callNextMethod()
})

#' Accessor for assignmentTable, not exported yet
#' @param object TopDownExperiment
#' @return data.table
#' @noRd
setMethod("assignmentTable", "TopDownExperiment", function(object) {
  object@assignmentTable
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
    prefix <- sprintf("Amino acid sequence (%d):", nchar(object@sequence))
    cat(prefix, .snippet(object@sequence, getOption("width") - nchar(prefix)),
        "\n")
  }

  if (nrow(fragmentTable(object))) {
    cat("Number of theoretical Fragments:", nrow(fragmentTable(object)), "\n")
    fragments <- sort(unique(.fragmentTypes(object)))
    prefix <- sprintf("Fragment types (%d):", length(fragments))
    cat(prefix, .snippet(paste0(fragments, collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (nrow(assignmentTable(object))) {
    cat("Number of Spectra-Fragment-Assignments:",
        nrow(assignmentTable(object)), "\n")
  }
  invisible(NULL)
})

