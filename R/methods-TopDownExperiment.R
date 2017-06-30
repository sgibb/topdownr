setMethod("show", "TopDownExperiment", function(object) {
  cat("TopDown Experiment data (\"", class(object), "\")\n", sep="")

  cat("dim:", dim(object), "\n")

  if (nchar(object@sequence)) {
    cat("Amino acid sequence:", object@sequence, "\n")
  }

  if (nrow(object@fragmentTable)) {
    cat("Number of theoretical Fragments:", nrow(object@fragmentTable), "\n")
    cat("Fragment types: ",
        paste0(sort(unique(object@fragmentTable$type)), collapse=", "), "\n")
  }

  if (nrow(object@assignmentTable)) {
    cat("Number of Spectra-Fragment-Assignments:",
        nrow(object@assignmentTable), "\n")
  }

  callNextMethod(object)
})
