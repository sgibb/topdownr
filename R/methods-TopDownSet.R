#' @param x TopDownSet
#' @noRd
setMethod("[", c("TopDownSet", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE) {
  d0 <- dim(x)
  dn <- dimnames(x)

  if (missing(i)) {
    i <- seq_len(d0[1L])
  }

  if (is.character(i)) {
    ii <- .subset(dn[[1L]] %in% i | as.character(fragmentType(x)) %in% i,
                  d0[1L], dn[[1L]])
  } else {
    ii <- .subset(i, d0[1L], dn[[1L]])
  }

  if (is.unsorted(ii)) {
    warning("It is not possible to change the row order.")
    ii <- sort(ii)
  }

  if (missing(j)) {
    j <- seq_len(d0[2L])
  }
  jj <- .subset(j, d0[2L], dn[[2L]])

  if (drop) {
    warning("'drop' is ignored.")
  }

  x@assay <- x@assay[ii, jj, drop=FALSE]
  x@colData <- .droplevels(x@colData[jj, ])
  x@rowViews <- x@rowViews[ii, ]
  x@rowViews@elementMetadata <- .droplevels(x@rowViews@elementMetadata)
  isFasta <- grepl(.topDownFileExtRx("fasta"), x@files)
  x@files <- x@files[.subsetFiles(x@files, unique(x@colData$File)) | isFasta]

  d1 <- dim(x)

  x <- .tdsLogMsg(x, paste0("Subsetted [", d0[1L], ";", d0[2L], "] to [",
                                           d1[1L], ";", d1[2L], "]."))

  if (validObject(x)) {
    x
  }
})

#' @param x TopDownSet
#' @noRd
setMethod("[[", c("TopDownSet", "ANY", "missing"), function(x, i, j, ...) {
  colData(x)[[i, ...]]
})

#' @param x TopDownSet
#' @noRd
setReplaceMethod("[[", c("TopDownSet", "ANY", "missing"),
                 function(x, i, j, ..., value) {
    colData(x)[[i, ...]] <- value
    x
})

#' @param x TopDownSet
#' @noRd
.DollarNames.TopDownSet <- function(x, pattern="")
    grep(pattern, names(colData(x)), value=TRUE)

#' @param x TopDownSet
#' @noRd
setMethod("$", "TopDownSet", function(x, name) {
    colData(x)[[name]]
})

#' @param x TopDownSet
#' @noRd
setReplaceMethod("$", "TopDownSet", function(x, name, value) {
    colData(x)[[name]] <- value
    x
})

#' @param object TopDownSet
#' @return numeric
#' @noRd
setMethod("aggregate", "TopDownSet", function(x,
                                              by=c("Mz", "AGCTarget",
                                                   "ETDReagentTarget",
                                                   "ETDActivation",
                                                   "CIDActivation",
                                                   "HCDActivation",
                                                   "SupplementalActivationCE",
                                                   "SupplementalActivation")) {
  d0 <- dim(x)

  if (!(all(by %in% colnames(x@colData)))) {
    stop(by[!by %in% colnames(x@colData)], " is/are not present in 'colData(x)'.")
  }

  groups <- .groupByLabels(as.data.frame(x@colData), by)
  x@assay <- .rowMeansGroup(x@assay, groups)
  x@colData <- .aggregateDataFrame(x@colData, groups,
                                   ignoreNumCols=c("Scan", "ConditionId"))
  ## now meaningless
  x@files <- x@files[grepl(.topDownFileExtRx("fasta"), x@files)]

  d1 <- dim(x)

  x <- .tdsLogMsg(x, paste0("Aggregated [", d0[1L], ";", d0[2L], "] to [",
                                            d1[1L], ";", d1[2L], "]."))

  if (validObject(x)) {
    x
  }
})

#' @param object TopDownSet
#' @return Matrix
#' @export
#' @noRd
setMethod("assayData", "TopDownSet", function(object) {
  object@assay
})

#' @param object TopDownSet
#' @return DataFrame
#' @export
#' @noRd
setMethod("colData", "TopDownSet", function(object) {
  object@colData
})

#' @param x TopDownSet
#' @return TopDownSet
#' @noRd
setReplaceMethod("colData", "TopDownSet", function(object, ..., value) {
  object@colData <- value
  object
})

#' @param object TopDownSet
#' @return DataFrame
#' @export
#' @noRd
setMethod("conditionData", "TopDownSet", function(object, ...) {
  colData(object)
})

#' @param object TopDownSet
#' @return TopDownSet
#' @export
#' @noRd
setReplaceMethod("conditionData", "TopDownSet", function(object, ..., value) {
  colData(object) <- value
  object
})

#' @param object TopDownSet
#' @return numeric
#' @noRd
setMethod("dim", "TopDownSet", function(x) {
  dim(x@assay)
})

#' @param object TopDownSet
#' @return list
#' @noRd
setMethod("dimnames", "TopDownSet", function(x) {
  list(names(x@rowViews), row.names(x@colData))
})

#' @param object TopDownSet
#' @return FragmentViews
#' @export
#' @noRd
setMethod("fragmentData", "TopDownSet", function(object, ...) {
  rowViews(object)
})

#' @param object TopDownSet
#' @return FragmentViews
#' @export
#' @noRd
setMethod("rowViews", "TopDownSet", function(object, ...) {
  object@rowViews
})

#' @param object TopDownExperiment
#' @noRd
setMethod("show", "TopDownSet", function(object) {
  cat(sprintf("%s object (%.2f Mb)\n",
              class(object), object.size(object) / 1024L^2L))

  if (length(object@rowViews)) {
    cat("- - - Protein data - - -\n")
    prefix <- sprintf("Amino acid sequence (%d):",
                      nchar(object@rowViews@subject))
    cat(prefix, .snippet(as.character(object@rowViews@subject),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (length(object@rowViews)) {
    cat("- - - Fragment data - - -\n")
    cat("Number of theoretical fragments:", length(object@rowViews), "\n")
    fragments <- fragmentType(object)
    cat0("Theoretical fragment types (", nlevels(fragments), "): ",
         paste0(.hft(levels(fragments), n=5), collapse=", "), "\n")
    mass <- range(fragmentMass(object))
    cat(sprintf("Theoretical mass range: [%.2f;%.2f]\n", mass[1L], mass[2L]))
  }

  if (nrow(object@colData)) {
    cat("- - - Condition data - - -\n")
    cat("Number of conditions:", nrow(object@colData), "\n")
    cat0("Condition variables (", ncol(object@colData), "): ",
         paste0(.hft(colnames(object@colData), n=2), collapse=", "), "\n")
  }

  if (all(dim(object))) {
    cat("- - - Intensity data - - -\n")
    cat(sprintf("Size of array: %dx%d (%.2f%% != 0)\n",
                nrow(object@assay), ncol(object@assay),
                length(object@assay@x) / length(object@assay) * 100L))
    if (length(object@assay@x)) {
      intensity <- range(object@assay@x)
    } else {
      intensity <- c(-NA, NA)
    }
    cat(sprintf("Intensity range: [%.2f;%.2f]\n", intensity[1L], intensity[2L]))
  }

  if (length(object@processing)) {
    cat("- - - Processing information - - -\n")
    cat(paste0(object@processing, collapse="\n"), "\n")
  }

  invisible(NULL)
})
