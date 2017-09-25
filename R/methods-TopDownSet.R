#' @describeIn TopDownSet Aggregate conditions/runs.
#'
#' Aggregates conditions/runs (columns) in an [TopDownSet-class] object `by` a
#' user-given value (default is the `"Sample"` column of `colData` which has the
#' same value for technical replicates).
#' It combines intensity values and numeric metadata of the grouped
#' conditions/runs (columns) by `mean` and returns a reduced
#' [TopDownSet-class] object.
#'
#' @param x,object `TopDownSet`
#' @param by `list`, grouping variable, in general it refers to technical
#' @param \ldots arguments passed to internal/other methods.
#' replicates (that's why the default is the `"Sample"` column in `colData`).
## @return `TopDownSet`
#' @aliases aggregate,TopDownSet-method
#' @export
setMethod("aggregate", "TopDownSet",
          function(x, by=x$Sample, ...) {
    d0 <- .logdim(x)
    groups <- .groupByLabels(by)

    if (length(groups) != ncol(x)) {
        stop("'by' has to be of the same length as 'ncol(x)'.")
    }
    x@assay <- .rowMeansGroup(x@assay, groups)
    x@colData <- .aggregateDataFrame(x@colData, groups,
                                     ignoreNumCols=c("Scan", "Condition"))

    d1 <- .logdim(x)

    x <- .atdsLogMsg(x, "Aggregated ", d0, " to ", d1, ".", addDim=FALSE)

    if (validObject(x)) {
        x
    }
})

#' @describeIn TopDownSet Filter by CV.
#'
#' Filtering is done by coefficient of variation across technical replicates
#' (defined by the `by` argument). All fragments below a given `threshold`
#' are removed. The `threshold` is the maximal allowed CV in percent (`sd/mean *
#' 100 < threshold`).
#'
## @param object `TopDownSet`
#' @param threshold `double`, threshold variable.
## @param by `list`, how technical replicates are defined.
## @return `TopDownSet`
#' @aliases filterCv filterCv,TopDownSet-method
#' @export
setMethod("filterCv", "TopDownSet",
          function(object, threshold, by=object$Sample, ...) {
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
        object <- .atdsLogMsg(object, n0 - n1, " fragments with CV > ",
                              threshold, "% filtered")

    }
    if (validObject(object)) {
        object
    }
})

#' @describeIn TopDownSet Filter by ion injection time.
#'
#' Filtering is done by maximal allowed deviation and just the technical
#' `keepTopN` replicates with the lowest deviation from the median ion
#' injection time are kept.
#'
## @param object `TopDownSet`
#' @param maxDeviation `double`, maximal allowed deviation in the `log2`
#' injection time in ms in comparison to the median ion injection time.
#' @param keepTopN `integer`, how many technical replicates should be kept?
## @param by `list`, how technical replicates are defined.
## @return `TopDownSet`
#' @aliases filterInjectionTime filterInjectionTime,TopDownSet-method
#' @export
setMethod("filterInjectionTime", "TopDownSet",
          function(object, maxDeviation=log2(3), keepTopN=2,
                   by=object$Sample, ...) {
    if (!is.numeric(maxDeviation) || !length(maxDeviation) == 1L) {
      stop("'maxDeviation' has to be a 'numeric' of length one.")
    }

    if (maxDeviation < 0) {
      stop("'maxDeviation' has to be greater than 0.")
    }

    lr <- as.vector(abs(log2(object$IonInjectionTimeMs /
                             object$MedianIonInjectionTimeMs)))

    i <- intersect(which(lr <= maxDeviation),
                   .topIdx(-lr, .groupByLabels(by), n=keepTopN))

    if (length(i) && length(i) != ncol(object)) {
        n0 <- ncol(object)
        object <- object[, i]
        n1 <- ncol(object)
        nd <- n0 - n1
        object <- .atdsLogMsg(object, n0 - n1, " scan", if (nd > 1L) { "s" },
                              " filtered with injection time deviation >= ",
                              maxDeviation, " or rank >= ", keepTopN + 1L)
    }
    if (validObject(object)) {
        object
    }
})

#' @describeIn TopDownSet Filter by intensity.
#'
#' Filtering is done by removing all fragments that are below a given
#' (absolute/relative) intensity `threshold`.
#
## @param object `TopDownSet`
## @param threshold `double`, remove fragments with intensity below
#' `threshold`.
#' @param relative `logical`, if relative is `TRUE` all fragments with
#' intensity below `threshold * max(intensity)` per fragment are removed,
#' otherwise all fragments below `threshold` are removed.
## @return `TopDownSet`
#' @aliases filterIntensity filterIntensity,TopDownSet-method
#' @export
setMethod("filterIntensity", "TopDownSet",
          function(object, threshold, relative=TRUE, ...) {
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
        object <- .atdsLogMsg(object, n0 - n1, " intensity values < ",
                              threshold, if (relative) { " (relative)" },
                              " filtered")
    }
    if (validObject(object)) {
        object
    }
})

#' @describeIn TopDownSet Filter by non-replicated fragments.
#'
#' Filtering is done by removing all fragments that don't replicate across
#' technical replicates.
#'
## @param object `TopDownSet`
#' @param minN `numeric`, if less than `minN` of a fragment are found across
#' technical replicates it is removed.
## @param by `list`, how technical replicates are defined.
## @return `TopDownSet`
#' @aliases filterNonReplicatedFragments
#' filterNonReplicatedFragments,TopDownSet-method
#' @export
setMethod("filterNonReplicatedFragments", "TopDownSet",
          function(object, minN=2, by=object$Sample, ...) {
    if (!is.numeric(minN) || !length(minN) == 1L) {
        stop("'minN' has to be a 'numeric' of length one.")
    }
    minN <- as.integer(minN)

    n0 <- nnzero(object@assay)

    object@assay <- .drop0rowReplicates(object@assay, .groupByLabels(by),
                                        minN=minN)
    n1 <- nnzero(object@assay)

    if (n0 - n1) {
        object <- .atdsLogMsg(object, n0 - n1, " intensity values of ",
                              "fragments replicated < ", minN,
                              " times filtered")
    }
    if (validObject(object)) {
        object
    }
})

#' @describeIn TopDownSet Normalise.
#'
#' Applies *T*otal *I*on *C*urrent normalisation to a [TopDownSet-class]. The
#' normalisation ist done per scans/conditions (column-wise normalisation).
#'
## @param object `TopDownSet`
#' @param method `character`, normalisation method, currently just `"TIC"`
#' for *T*otal *I*on *C*urrent normalisation of the scans/conditions
#' (column-wise normalisation) is supported.
## @return `TopDownSet`
#' @aliases normalize,TopDownSet-method
#' @export
setMethod("normalize", "TopDownSet",
          function(object, method="TIC", ...) {
    method <- match.arg(method)

    if (method == "TIC") {
        object@assay <- .normaliseCols(object@assay,
                                       scale=object$TotIonCurrent)
        .atdsLogMsg(object, "Intensity values normalized to TIC.",
                    addDim=FALSE)
    }
})

#' @rdname TopDownSet-class
#' @aliases show,TopDownSet-method
#' @export
setMethod("show", "TopDownSet", function(object) {
    callNextMethod()

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
        cat("Number of matched fragments:", nnzero(object@assay), "\n")
        cat(sprintf("Intensity range: [%.2f;%.2f]\n",
                    intensity[1L], intensity[2L]))
    }

    if (length(object@processing)) {
        cat("- - - Processing information - - -\n")
        cat(paste0(object@processing, collapse="\n"), "\n")
    }

    invisible(NULL)
})

#' @describeIn TopDownSet Summary statistics.
#'
#' Returns a `matrix` with some statistics: number of fragments,
#' total/min/first quartile/median/mean/third quartile/maximum of intensity
#' values.
#'
## @param object `TopDownSet`
#' @param what `character`, specifies whether `"conditions"` (columns; default) or
#' `"fragments"` (rows) should be summarized.
#' @aliases summary,TopDownSet-method
#' @export
setMethod("summary", "TopDownSet",
          function(object, what=c("conditions", "fragments"), ...) {
    what <- if (match.arg(what) == "conditions") { "columns" } else { "rows" }
    callNextMethod(object=object, what=what)
})

#' @rdname TopDownSet-class
#' @name as
#' @section Coercion:
#'
#' `as(object, "NCBSet"): Coerce an [TopDownSet-class] object into an
#' [NCBSet-class] object.
setAs("TopDownSet", "NCBSet", function(from) {
    assay <- .ncbMap(from)
    ncb <- new("NCBSet",
               rowViews=Views(subject(from@rowViews),
                              start=1L, width=seq_len(nrow(assay)),
                              names=rownames(assay)),
               colData=from@colData,
               assay=assay,
               files=from@files,
               tolerance=from@tolerance,
               processing=from@processing)
    .atdsLogMsg(ncb, "Coerced TopDownSet into an NCBSet object")
})
