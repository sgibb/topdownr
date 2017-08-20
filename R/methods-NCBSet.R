#' Find best combination of conditions for highest coverage
#'
#' @param object `NCBSet`
#' @param n `integer`, max number of combinations/iterations.
#' @param minN `integer` stop if there are less than `minN` additional fragments
#' added.
#' @return `matrix`, first column: index of condition, second column: number of
#' added fragments.
#' @export
#' @noRd
setMethod("bestConditions", "NCBSet",
          function(object, n=ncol(object), minN=0L) {
    m <- .bestCoverageCombination(object@assay, n=n, minN=minN)
    colnames(m) <- c("index", "fragments")
    m
})

#' @param object `NCBSet`
#' @noRd
setMethod("show", "NCBSet", function(object) {
    callNextMethod()

    if (length(object@rowViews)) {
        cat("- - - Fragment data - - -\n")
        cat0(paste0("Number of ",
                    c("N-terminal", "C-terminal", "N- and C-terminal"),
                    " fragments: ", tabulate(object@assay@x, nbins=3L), "\n"))
    }

    if (nrow(object@colData)) {
        cat("- - - Condition data - - -\n")
        cat("Number of conditions:", length(unique(object$Sample)), "\n")
        cat("Number of scans:", nrow(object@colData), "\n")
        cat0("Condition variables (", ncol(object@colData), "): ",
             paste0(.hft(colnames(object@colData), n=2), collapse=", "), "\n")
    }

    if (all(dim(object))) {
        cat("- - - Assay data - - -\n")
        cat(sprintf("Size of array: %dx%d (%.2f%% != 0)\n",
                    nrow(object@assay), ncol(object@assay),
                    nnzero(object@assay) / length(object@assay) * 100L))
    }

    if (length(object@processing)) {
        cat("- - - Processing information - - -\n")
        cat(paste0(object@processing, collapse="\n"), "\n")
    }

    invisible(NULL)
})
