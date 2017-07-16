#' @param x TopDownSet
#' @noRd
setMethod("[", "TopDownSet", function(x, i, j="missing", ...,
                                             drop="missing") {
  stop("NOT IMPLEMENTED YET")
})

#' @param object TopDownSet
#' @return numeric
#' @noRd
setMethod("dim", "TopDownSet", function(x) {
  dim(x@assays)
})

#' @param object TopDownSet
#' @return list
#' @noRd
setMethod("dimnames", "TopDownSet", function(x) {
  list(names(x@rowViews), row.names(x@colData))
})

#' @param object TopDownExperiment
#' @noRd
setMethod("show", "TopDownSet", function(object) {
  cat(sprintf("%s object (%.2f Mb)\n",
              class(object), object.size(object) / 1024L^2L))

  if (nchar(object@rowViews@subject)) {
    cat("- - - Protein data - - -\n")
    prefix <- sprintf("Amino acid sequence (%d):",
                      nchar(object@rowViews@subject))
    cat(prefix, .snippet(object@rowViews@subject,
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (length(object@rowViews)) {
    cat("- - - Fragment data - - -\n")
    cat("Number of theoretical fragments:", length(object@rowViews), "\n")
    fragments <- fragmentType(object)
    prefix <- sprintf("Theoretical fragment types (%d):", nlevels(fragments))
    cat(prefix, .snippet(paste0(levels(fragments), collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
    mass <- range(fragmentMass(object))
    cat(sprintf("Theoretical mass range: [%.2f;%.2f]\n",
                mass[1L], mass[2L]))
  }

  if (nrow(object@colData)) {
    cat("- - - Condition data - - -\n")
    cat("Number of conditions:", nrow(object@colData), "\n")
    prefix <- sprintf("Condition variables (%d):", ncol(object@colData))
    cat(prefix, .snippet(paste0(colnames(object@colData), collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (all(dim(object))) {
    cat("- - - Intensity data - - -\n")
    cat(sprintf("Size of array: %dx%d (%.2f%% != 0)\n", nrow(object@assays),
        ncol(object@assays), mean(object@assays != 0L) * 100L))
    intensity <- range(object@assays@x)
    cat(sprintf("Intensity range: [%.2f;%.2f]\n", intensity[1L], intensity[2L]))
  }

  if (length(object@processing)) {
    cat("- - - Processing information - - -\n")
    cat(paste0(object@processing, collapse="\n"), "\n")
  }

  invisible(NULL)
})

