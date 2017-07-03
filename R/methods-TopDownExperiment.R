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
  cat("TopDown Experiment data (\"", class(object), "\")\n", sep="")

  if (nchar(object@sequence)) {
    cat("Amino acid sequence:", object@sequence, "\n")
  }

  if (nrow(fragmentTable(object))) {
    cat("Number of theoretical Fragments:", nrow(fragmentTable(object)), "\n")
    cat("Fragment types: ",
        paste0(sort(unique(.fragmentTypes(object))), collapse=", "), "\n")
  }

  if (nrow(assignmentTable(object))) {
    cat("Number of Spectra-Fragment-Assignments:",
        nrow(assignmentTable(object)), "\n")
  }
  invisible(NULL)
})

