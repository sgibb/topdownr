#' @param x [TopDownSet-class]
#' @param i `numeric`, `logical` or `character`, subsetting on fragment data,
#' names (`c("a1", "b1", "c1", "c2", "c3")`) and types (`c("c", "x")`) are
#' supported.
#' @param j `numeric` or `logical`, subsetting based on condition data.
#' @param \ldots currently ignored.
#' @param drop `logical`, currently ignored.
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

   x <- .tdsLogMsg(x, "Subsetted [", d0[1L], ";", d0[2L], "] to [",
                                     d1[1L], ";", d1[2L], "].")

   if (validObject(x)) {
       x
   }
})

#' @param x [TopDownSet-class]
#' @param i `numeric`, `logical` or `character`, subsetting based on condition
#' data.
#' @param j currently ignored.
#' @param \ldots currently ignored.
#' @noRd
setMethod("[[", c("TopDownSet", "ANY", "missing"), function(x, i, j, ...) {
    colData(x)[[i, ...]]
})

#' @param x `TopDownSet`
#' @noRd
setReplaceMethod("[[", c("TopDownSet", "ANY", "missing"),
                 function(x, i, j, ..., value) {
    colData(x)[[i, ...]] <- value
    if (validObject(x)) {
      x
    }
})

#' @param x `TopDownSet`
#' @noRd
#' @export
.DollarNames.TopDownSet <- function(x, pattern="") {
    grep(pattern, names(colData(x)), value=TRUE)
}

#' @param x `TopDownSet`
#' @noRd
setMethod("$", "TopDownSet", function(x, name) {
    colData(x)[[name]]
})

#' @param x `TopDownSet`
#' @noRd
setReplaceMethod("$", "TopDownSet", function(x, name, value) {
    colData(x)[[name]] <- value
    x
})

#' @param object `TopDownSet`
#' @param by `list`, grouping variable
#' @return `TopDownSet`
#' @noRd
setMethod("aggregate", "TopDownSet",
          function(x, by=x$Sample) {
    d0 <- dim(x)
    groups <- .groupByLabels(by)

    if (length(groups) != ncol(x)) {
        stop("'by' has to be of the same length as 'ncol(x)'.")
    }
    x@assay <- .rowMeansGroup(x@assay, groups)
    x@colData <- .aggregateDataFrame(x@colData, groups,
                                     ignoreNumCols=c("Scan", "Condition"))
    ## now meaningless
    x@files <- x@files[grepl(.topDownFileExtRx("fasta"), x@files)]

    d1 <- dim(x)

    x <- .tdsLogMsg(x, "Aggregated [", d0[1L], ";", d0[2L], "] to [",
                                       d1[1L], ";", d1[2L], "].")

    if (validObject(x)) {
        x
    }
})

#' @param object `TopDownSet`
#' @return `Matrix`
#' @export
#' @noRd
setMethod("assayData", "TopDownSet", function(object) {
    object@assay
})

#' @param object `TopDownSet`
#' @return `DataFrame`
#' @export
#' @noRd
setMethod("colData", "TopDownSet", function(object) {
    object@colData
})

#' @param x `TopDownSet`
#' @return `TopDownSet`
#' @export
#' @noRd
setReplaceMethod("colData", "TopDownSet", function(object, ..., value) {
    object@colData <- value
    if (validObject(object)) {
      object
    }
})

#' @param object `TopDownSet`
#' @return `DataFrame`
#' @export
#' @noRd
setMethod("conditionData", "TopDownSet", function(object, ...) {
    colData(object)
})

#' @param object `TopDownSet`
#' @return `TopDownSet`
#' @export
#' @noRd
setReplaceMethod("conditionData", "TopDownSet", function(object, ..., value) {
    colData(object) <- value
    if (validObject(object)) {
      object
    }
})

#' @param object `TopDownSet`
#' @return `numeric`
#' @noRd
setMethod("dim", "TopDownSet", function(x) {
    dim(x@assay)
})

#' @param object `TopDownSet`
#' @return `list`
#' @noRd
setMethod("dimnames", "TopDownSet", function(x) {
    list(names(x@rowViews), row.names(x@colData))
})

#' @param object `TopDownSet`
#' @return `FragmentViews`
#' @export
#' @noRd
setMethod("fragmentData", "TopDownSet", function(object, ...) {
    rowViews(object)
})

#' Filter `TopDownSet` by CV.
#'
#' Filtering is done by coefficent of variation across technical replicates.
#' All fragments below a given threshold are removed.
#'
#' @param object `TopDownSet`
#' @param threshold max allowed CV in percent (`sd/mean * 100 < threshold`).
#' @param by `list`, how technical repliactes are defined.
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("filterCv", "TopDownSet",
          function(object, threshold, by=object$Sample) {
    if (!is.numeric(threshold) || !length(threshold) == 1L) {
      stop("'threshold' has to be a 'numeric' of length one.")
    }

    if (threshold < 0) {
      stop("'threshold' has to be greater than 0.")
    }

    n0 <- nnzero(object@assay)

    group <- .groupByLabels(by)
    mm <- .createMaskMatrix(group)
    cvs <- tcrossprod(.rowCvsGroup(object@assay, group, na.rm=TRUE), mm)
    object@assay[Matrix::which(cvs * 100L > threshold)] <- 0L
    object@assay <- drop0(object@assay, is.Csparse=TRUE)

    n1 <- nnzero(object@assay)
    if (n0 - n1) {
        object <- .tdsLogMsg(object, n0 - n1, " fragments with CV > ",
                             threshold, "% filtered.")

    }
    if (validObject(object)) {
        object
    }
})

#' Filter `TopDownSet` by ion injection time.
#'
#' Filtering is done by maximal allowed deviation and just the technical
#' `keepTopN` replicates with the lowest deviation from the median ion
#' injection time are kept.
#'
#' @param object `TopDownSet`
#' @param maxDeviation `double`, maximal allowed deviation in the log2 injection
#' time in ms in comparison to the median ion injection time.
#' @param keepTopN `integer`, how many technical repliactes should be kept?
#' @param by `list`, how technical repliactes are defined.
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("filterInjectionTime", "TopDownSet",
          function(object, maxDeviation=log2(3), keepTopN=2,
                   by=object$Sample) {
    if (!is.numeric(maxDeviation) || !length(maxDeviation) == 1L) {
      stop("'maxDeviation' has to be a 'numeric' of length one.")
    }

    if (maxDeviation < 0) {
      stop("'maxDeviation' has to be greater than 0.")
    }

    lr <- log2(object$IonInjectionTimeMs / object$MedianIonInjectionTimeMs)

    i <- intersect(which(lr <= maxDeviation),
                   .topIdx(-abs(as.vector(lr)), .groupByLabels(by),
                           n=keepTopN))

    if (length(i) && length(i) != ncol(object)) {
        n0 <- ncol(object)
        object <- object[, i]
        n1 <- ncol(object)
        nd <- n0 - n1
        object <- .tdsLogMsg(object, n0 - n1, " scan", if (nd > 1L) { "s" },
                             " filtered with injection time deviation >= ",
                             maxDeviation, " or rank >= ", keepTopN + 1L, ".")
    }
    if (validObject(object)) {
        object
    }
})

#' Filter `TopDownSet` by intensity.
#'
#' Filtering is done by removing all fragments that are below a given
#' (absolute/relative) intensity threshold.
#
#' @param object `TopDownSet`
#' @param threshold `double`, remove fragments with intensity below
#' `threshold`.
#' @param relative `logical`, if relative is `TRUE` all fragments with
#' intensity below `threshold * max(intensity)` per fragment are removed.
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("filterIntensity", "TopDownSet",
          function(object, threshold, relative=TRUE) {
    if (!is.numeric(threshold) || !length(threshold) == 1L) {
        stop("'threshold' has to be a 'numeric' of length one.")
    }

    n0 <- nnzero(object@assay)

    if (relative) {
        if (1L < threshold || threshold < 0L) {
            stop("'threshold hast to be between 0 and 1.")
        }
        object@assay <- .drop0rowLt(object@assay,
                                    tol=.rowMax(object@assay) * threshold)
    } else {
        object@assay <- drop0(object@assay,
                              tol=threshold - 10L * .Machine$double.eps,
                              is.Csparse=TRUE)
    }
    n1 <- nnzero(object@assay)
    if (n0 - n1) {
        object <- .tdsLogMsg(object, n0 - n1, " intensity values < ",
                             threshold, if (relative) { " (relative)" },
                             " filtered.")
    }
    if (validObject(object)) {
        object
    }
})

#' Filter `TopDownSet` by intensity.
#'
#' Filtering is done by removing all fragments that don't replicate across
#' technical replicates.
#'
#' @param object `TopDownSet`
#' @param minN `numeric`, if less than `minN` of a fragment are found across
#' technical replicates it is removed.
#' @param by `list`, how technical repliactes are defined.
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("filterNonReplicatedFragments", "TopDownSet",
          function(object, minN=2, by=object$Sample) {
    if (!is.numeric(minN) || !length(minN) == 1L) {
        stop("'minN' has to be a 'numeric' of length one.")
    }
    minN <- as.integer(minN)

    n0 <- nnzero(object@assay)

    object@assay <- .drop0rowReplicates(object@assay, .groupByLabels(by),
                                        minN=minN)
    n1 <- nnzero(object@assay)

    if (n0 - n1) {
        object <- .tdsLogMsg(object, n0 - n1, " intensity values of ",
                             "fragments replicated < ", minN,
                             " times filtered.")
    }
    if (validObject(object)) {
        object
    }
})

#' @param object `TopDownSet`
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("normalize", "TopDownSet",
          function(object, ...) {
    object@assay <- .normaliseRows(object@assay)
    .tdsLogMsg(object, "Fragment intensity values normalized.")
})

#' @param object `TopDownSet`
#' @return `TopDownSet`
#' @export
#' @noRd
setMethod("removeEmptyConditions", "TopDownSet",
          function(object) {
    i <- Matrix::colSums(object@assay) != 0L
    object <- object[, i]
    .tdsLogMsg(object, sum(!i), " empty conditions removed.")
})

#' @param object `TopDownSet`
#' @return `FragmentViews`
#' @export
#' @noRd
setMethod("rowViews", "TopDownSet", function(object, ...) {
    object@rowViews
})

#' @param object `TopDownSet`
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
        if (length(metadata(object@rowViews)$modifications)) {
            mod <- metadata(object@rowViews)$modifications
            cat0("Modifications (", length(mod), "): ",
                 paste0(.hft(mod, n=3), collapse=", "), "\n")
        }
    }

    if (length(object@rowViews)) {
        cat("- - - Fragment data - - -\n")
        cat("Number of theoretical fragments:", length(object@rowViews), "\n")
        fragments <- fragmentType(object)
        cat0("Theoretical fragment types (", nlevels(fragments), "): ",
             paste0(.hft(levels(fragments), n=5), collapse=", "), "\n")
        mass <- range(fragmentMass(object))
        cat(sprintf("Theoretical mass range: [%.2f;%.2f]\n",
                    mass[1L], mass[2L]))
    }

    if (nrow(object@colData)) {
        cat("- - - Condition data - - -\n")
        cat("Number of conditions:", length(unique(object$Sample)), "\n")
        cat("Number of scans:", nrow(object@colData), "\n")
        cat0("Condition variables (", ncol(object@colData), "): ",
             paste0(.hft(colnames(object@colData), n=2), collapse=", "), "\n")
    }

    if (all(dim(object))) {
        cat("- - - Intensity data - - -\n")
        cat(sprintf("Size of array: %dx%d (%.2f%% != 0)\n",
                    nrow(object@assay), ncol(object@assay),
                    nnzero(object@assay) / length(object@assay) * 100L))
        if (length(object@assay@x)) {
            intensity <- range(object@assay@x)
        } else {
            intensity <- c(-NA, NA)
        }
        cat(sprintf("Intensity range: [%.2f;%.2f]\n",
                    intensity[1L], intensity[2L]))
    }

    if (length(object@processing)) {
        cat("- - - Processing information - - -\n")
        cat(paste0(object@processing, collapse="\n"), "\n")
    }

    invisible(NULL)
})
