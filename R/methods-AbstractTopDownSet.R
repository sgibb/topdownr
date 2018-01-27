#' @describeIn AbstractTopDownSet Subset operator.
#'
#' For `i` `numeric`, `logical` or `character` vectors or empty
#' (missing) or `NULL` are supported.
#' Subsetting is done on the fragment/bond (row) level.
#' `character` indices could be names
#' (e.g. `c("a1", "b1", "c1", "c2", "c3")`)
#' or types (e.g. `c("c", "x")`) of the fragments for
#' [TopDownSet-class] objects,
#' or names of the bonds (e.g. `c("bond001")`) for
#' [NCBSet-class] objects. \cr
#' `j` could be a `numeric` or `logical` vector
#' and subsetting is done on the condition/run (column) level.
#'
#' @param object,x `AbstractTopDownSet`
#' @param i,j `numeric`, `logical` or `character`,
#' indices specifying elements to extract or replace.
#' @param \ldots arguments passed to internal/other methods.
#' @param drop `logical`, currently ignored.
#'
## @param i `numeric`, `logical` or `character`, subsetting on fragment/bond
## data, names (`c("a1", "b1", "c1", "c2", "c3")`) and types (`c("c", "x")`)
## (or for [NCBSet-class] bonds e.g. `bond001`) are supported.
## @param j `numeric` or `logical`, subsetting based on condition data.
## @param \ldots currently ignored.'
#' @aliases [,AbstractTopDownSet,ANY,ANY,ANY-method
#' @export
setMethod("[", c("AbstractTopDownSet", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE) {
    d0 <- dim(x)
    ld0 <- .logdim(x)
    dn <- dimnames(x)

    if (missing(i)) {
        i <- seq_len(d0[1L])
    }

    if (is.character(i)) {
        sel <- dn[[1L]] %in% i
        fragtype <- fragmentType(x)
        if (!is.null(fragtype)) {
            sel <- sel | as.character(fragtype) %in% i
        }
        ii <- .subset(sel, d0[1L], dn[[1L]])
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

    if (!is.null(x@rowViews@elementMetadata)) {
        x@rowViews@elementMetadata <- .droplevels(x@rowViews@elementMetadata)
    }
    isFasta <- grepl(.topDownFileExtRx("fasta"), x@files)
    x@files <- x@files[
        .subsetFiles(basename(x@files), unique(x@colData$File)) | isFasta
    ]

    ld1 <- .logdim(x)

    if (validObject(x)) {
        .atdsLogMsg(x, "Subsetted ", ld0, " to ", ld1, ".", addDim=FALSE)
    }
})

#' @describeIn AbstractTopDownSet Subset operator.
#'
#' `i` could be a `numeric` or `logical` vector and
#' subsetting is done on the condition/run (column) level.
#'
## @param x `AbstractTopDownSet`
## @param i `logical` or `character`, subsetting based on condition
## @param j currently ignored.
## @param \ldots currently ignored.
#' @aliases [[,AbstractTopDownSet,ANY,missing,-method
#' @export
setMethod("[[", c("AbstractTopDownSet", "ANY", "missing"),
          function(x, i, j, ...) {
    colData(x)[[i, ...]]
})

#' @describeIn AbstractTopDownSet Setter for a column in the `colData` slot.
#'
#' The `[[<-` operator is used to add/replace
#' a single column of the `colData` `DataFrame`.
#
## @param x `AbstractTopDownSet`
## @param i `logical` or `character`, subsetting based on condition
## @param j currently ignored.
## @param \ldots currently ignored.
#' @param value replacment value.
#' @aliases [[<-,AbstractTopDownSet,ANY,missing,-method
#' @export
setReplaceMethod("[[", c("AbstractTopDownSet", "ANY", "missing"),
                 function(x, i, j, ..., value) {
    colData(x)[[i, ...]] <- value
    if (validObject(x)) {
        x
    }
})

#' @param x `AbstractTopDownSet`
#' @export
#' @noRd
.DollarNames.AbstractTopDownSet <- function(x, pattern="") {
    grep(pattern, names(colData(x)), value=TRUE)
}

#' @describeIn AbstractTopDownSet Accessor for columns in the `colData` slot.
#'
#' The `$` simplifies the accession of a single column
#' of the `colData`.
#' It is identical to the `[[` operator.
#'
## @param x `AbstractTopDownSet`
#' @param name `character` name of an (non)existing column in `colData`.
## @return `AbstractTopDownSet`
#' @aliases $,AbstractTopDownSet-method
#' @export
setMethod("$", "AbstractTopDownSet", function(x, name) {
    colData(x)[[name]]
})

#' @describeIn AbstractTopDownSet Setter for a column in the `colData` slot.
#'
#' The `$<-` operator is used to add/replace a single column
#' of the `colData` `DataFrame`.
#' It is identical to the `[[<-` operator.
#'
## @param x `AbstractTopDownSet`
## @return `AbstractTopDownSet`
#' @aliases $<-,AbstractTopDownSet-method
#' @export
setReplaceMethod("$", "AbstractTopDownSet", function(x, name, value) {
    colData(x)[[name]] <- value
    x
})

#' @describeIn AbstractTopDownSet Test if two objects are equal.
#'
#' @param target, current `AbstractTopDownSet`
#' @return Either 'TRUE' or a the return value of `all.equal.default`.
#' @noRd
#' @export
all.equal.AbstractTopDownSet <- function(target, current, ..., check.date=FALSE) {
    if (!check.date) {
        target@processing <- gsub("^\\[[^]]+\\] *", "", target@processing)
        current@processing <- gsub("^\\[[^]]+\\] *", "", current@processing)
    }
    all.equal.default(target, current, ...)
}

#' @describeIn AbstractTopDownSet Accessor for the `assay` slot.
#'
#' Returns a [Matrix::dgCMatrix-class] that stores the
#' intensity/coverage information of [AbstractTopDownSet-class]
#' object.
#'
## @param object `AbstractTopDownSet`
## @return `dgCMatrix`
#' @aliases assayData,AbstractTopDownSet-method
#' @export
setMethod("assayData", "AbstractTopDownSet", function(object) {
    object@assay
})

#' @describeIn AbstractTopDownSet Accessor for the `colData` slot.
#'
#' Returns a [S4Vectors::DataFrame-class] that stores
#' metadata for the conditons/runs (columns) of the
#' [AbstractTopDownSet-class] object.
#'
## @param object `AbstractTopDownSet`
## @return `DataFrame`
#' @aliases colData colData,AbstractTopDownSet-method
#' @export
setMethod("colData", "AbstractTopDownSet", function(object) {
    object@colData
})

#' @describeIn AbstractTopDownSet Setter for the `colData` slot.
#'
#' Replaces metadata for the conditons/runs (columns) of the
#' [AbstractTopDownSet-class] object.
#'
## @param object `AbstractTopDownSet`
## @return `AbstractTopDownSet`
#' @aliases colData<- colData<-,AbstractTopDownSet-method
#' @export
setReplaceMethod("colData", "AbstractTopDownSet", function(object, ..., value) {
    object@colData <- value
    if (validObject(object)) {
        object
    }
})

#' @describeIn AbstractTopDownSet Combine `AbstractTopDownSet` objects.
#'
#' If the `rowViews` are identically `combine` allows to combine two or more
#' `AbstractTopDownSet` objects. Please note that it uses the default
#' `sampleColumns` to define technical replicates (see [readTopDownFiles()]).and
#' the default `by` argument to group ion injection times for the calculation of
#' the median time (see [updateMedianInjectionTime()]). Both could be modified
#' after `combine` by calling [updateConditionNames()] (with modified
#' `sampleColumns` argument) and [updateMedianInjectionTime()] (with modified
#' `by` argument).
#'
## @param x `AbstractTopDownSet`
#' @param y `AbstractTopDownSet`
## @param \ldots
#' @aliases combine combine,AbstractTopDownSet,AbstractTopDownSet-method
#' @export
setMethod("combine",
          signature(x="AbstractTopDownSet", y="AbstractTopDownSet"),
          function(x, y) {
    if (class(x) != class(y)) {
        stop(paste0("Objects must be the same class, but are ",
                    class(x), ", ", class(y), "."))
    }
    ldx0 <- .logdim(x)
    ldy0 <- .logdim(y)
    x@rowViews <- combine(x@rowViews, y@rowViews)
    x@colData <- .colsToRle(.colsToLogical(.rbind(x@colData, y@colData)))
    x@assay <- .cbind(x@assay, y@assay)[names(x@rowViews),]
    x@files <- unique(x@files, y@files)
    x@tolerance <- max(x@tolerance, y@tolerance)
    x@processing <- c(x@processing, y@processing)
    x <- updateConditionNames(x)
    x <- updateMedianInjectionTime(x)

    if (validObject(x)) {
        x <- .atdsLogMsg(
            x, "Combined ", ldx0, " and ", ldy0, " into a ", .logdim(x), " ",
            class(x), " object.", addDim=FALSE
        )
    }
})

#' @describeIn AbstractTopDownSet Accessor for the `colData` slot.
#'
#' An alias for `colData`.
#'
## @param object `AbstractTopDownSet`
## @return `DataFrame`
#' @aliases conditionData conditionData,AbstractTopDownSet-method
#' @export
setMethod("conditionData", "AbstractTopDownSet", function(object, ...) {
    colData(object)
})

#' @describeIn AbstractTopDownSet Setter for the `colData` slot.
#'
#' An alias for `colData<-`.
#'
## @param object `AbstractTopDownSet`
## @return `AbstractTopDownSet`
#' @aliases conditionData<- conditionData<-,AbstractTopDownSet-method
#' @export
setReplaceMethod("conditionData", "AbstractTopDownSet",
                 function(object, ..., value) {
    colData(object) <- value
    if (validObject(object)) {
        object
    }
})

#' @describeIn AbstractTopDownSet Accessor for condition names.
#'
#' Returns a `character` with names for the conditions/runs (columns).
#'
## @param object `AbstractTopDownSet`
## @return `character`
#' @aliases conditionNames conditionNames,AbstractTopDownSet-method
#' @export
setMethod("conditionNames", "AbstractTopDownSet", function(object) {
    dimnames(object)[[2L]]
})

#' @describeIn AbstractTopDownSet Accessor for dimensions.
#'
#' Returns a `numeric` with number of fragments/bonds (rows) and
#' conditions/runs (columns).
#'
## @param object `AbstractTopDownSet`
## @return `numeric`
#' @aliases dim,AbstractTopDownSet-method
#' @export
setMethod("dim", "AbstractTopDownSet", function(x) {
    dim(x@assay)
})

#' @describeIn AbstractTopDownSet Accessor for dimension names.
#'
#' Returns a `list` with names for the fragments/bonds (rows) and for the
#' conditions/runs (columns).
#'
## @param object `AbstractTopDownSet`
## @return `list`
#' @aliases dimnames,AbstractTopDownSet-method
#' @export
setMethod("dimnames", "AbstractTopDownSet", function(x) {
    list(names(x@rowViews), row.names(x@colData))
})

#' @describeIn AbstractTopDownSet Remove empty conditions/runs.
#'
#' Removes conditions/runs (columns) without any intensity/coverage
#' information from the [AbstractTopDownSet-class] object.
#' It returns a modified [AbstractTopDownSet-class] object.
#'
## @param object `AbstractTopDownSet`
## @return `AbstractTopDownSet`
#' @aliases removeEmptyConditions
#' removeEmptyConditions,AbstractTopDownSet-method
#' @export
setMethod("removeEmptyConditions", "AbstractTopDownSet",
          function(object) {
    i <- Matrix::colSums(object@assay) != 0L
    object <- object[, i]
    .atdsLogMsg(object, sum(!i), " empty conditions removed")
})

#' @describeIn AbstractTopDownSet Accessor for the `rowViews` slot.
#'
#' Depending on the implementation it returns an
#' [FragmentViews-class] object for
#' [TopDownSet-class] objects
#' or an [Biostrings::XStringViews] object for
#' [NCBSet-class] objects.
#'
## @param object `AbstractTopDownSet`
## @return `XStringViews`
#' @aliases rowViews rowViews,AbstractTopDownSet-method
#' @export
setMethod("rowViews", "AbstractTopDownSet", function(object, ...) {
    object@rowViews
})


#' @rdname AbstractTopDownSet-class
#' @aliases show,AbstractTopDownSet-method
#' @export
setMethod("show", "AbstractTopDownSet", function(object) {
    cat(sprintf("%s object (%.2f Mb)\n",
                class(object), object.size(object) / 1024L^2L))

    if (length(object@rowViews)) {
        cat("- - - Protein data - - -\n")
        prefix <- sprintf(
            "Amino acid sequence (%d):", nchar(object@rowViews@subject)
        )
        cat(prefix,
            .snippet(
                as.character(object@rowViews@subject),
                getOption("width") - nchar(prefix)
            ), "\n"
        )
        if (length(metadata(object@rowViews)$mass)) {
            cat0("Mass : ", metadata(object@rowViews)$mass, "\n")
        }
        if (length(metadata(object@rowViews)$modifications)) {
            mod <- metadata(object@rowViews)$modifications
            cat0("Modifications (", length(mod), "): ",
                paste0(.hft(mod, n=3), collapse=", "), "\n")
        }
    }

    invisible(NULL)
})

#' @describeIn AbstractTopDownSet Summary statistics.
#'
#' Returns a `matrix` with some statistics: number of fragments,
#' total/min/first quartile/median/mean/third quartile/maximum of intensity
#' values.
#'
## @param object `AbstractTopDownSet`
#' @param what `character`,
#' specifies whether `"rows"` or `"columns"` should be summarized.
#' @aliases summary,AbstractTopDownSet-method
#' @export
setMethod("summary", "AbstractTopDownSet",
          function(object, what=c("rows", "columns"), ...) {
    .summary(object@assay, what=what)
})

#' @describeIn AbstractTopDownSet Update condition names.
#'
#' Updates condition names based on `sampleColumns` from
#' `conditionData`/`colData`. Columns with just identical entries are ignored.
#' This method will create/update the `colData(object)$Sample` column that
#' identifies technical replicates and could be used in other methods.
#'
## @param object `AbstractTopDownSet`
#' @param sampleColumns `character`,
#' column names of the [colData()]
#' used to define a sample (technical replicate). This is used to add the
#' `Sample` column (used for easier aggregation, etc.).
#' @param verbose `logical`, verbose output?
#' @aliases updateConditionNames updateConditionNames,AbstractTopDownSet-method
#' @export
setMethod("updateConditionNames", "AbstractTopDownSet",
          function(object, sampleColumns=c("Mz", "AgcTarget",
                                           "EtdReagentTarget",
                                           "EtdActivation",
                                           "CidActivation",
                                           "HcdActivation",
                                           "UvpdActivation"),
                   verbose=interactive()) {

    isColPresent <- sampleColumns %in% colnames(object@colData)

    .msg(
        verbose & any(!isColPresent),
        paste0(sampleColumns[!isColPresent], collapse=", "),
        " is/are not present and ignored."
    )

    if (all(!isColPresent)) {
        stop("'sampleColumns' must contain names from 'colnames(colData(x))'")
    }

    sampleColumns <- sampleColumns[isColPresent]

    o <- .orderByColumns(object@colData, cols=sampleColumns)

    resorted <- is.unsorted(o)
    .msg(verbose & resorted, "Order of conditions changed.")

    object@colData <- object@colData[o, ]
    object@colData$Sample <- .groupId(object@colData, cols=sampleColumns)
    object@colData <- .colsToRle(object@colData)
    object@assay <- object@assay[, o]
    rn <- rownames(object@colData)
    colnames(object@assay) <-
        rownames(object@colData) <-
            .makeRowNames(object@colData[, sampleColumns, drop=FALSE])

    if (is.null(rn) || any(rn != rownames(object@colData))) {
        object <- .atdsLogMsg(
            object,
            "Condition names updated based on: ",
            paste0(sampleColumns, collapse=", "), ".",
            if (resorted) { " Order of conditions changed." },
            " ", length(unique(object@colData$Sample)), " conditions.",
            addDim=FALSE
        )
    }

    if (validObject(object)) {
        object
    }
})

#' @describeIn AbstractTopDownSet Update median ion injection times.
#'
#' Recalculates median ion injection times by a user given grouping variable
#' (default: Mz, AgcTarget). This is useful if you acquire new data and the ion
#' injection time differs across the runs. Use the `by` argument to provide a
#' `list`/`data.frame` of grouping variables, e.g.
#' `by=colData(object)[, c("Mz", "AgcTarget", "File")]`.
#'
## @param object `TopDownSet`
#' @param by `list`, grouping information.
## @return `AbstractTopDownSet`
#' @aliases updateMedianInjectionTime
#' updateMedianInjectionTime,TopDownSet-method
#' @export
setMethod("updateMedianInjectionTime", "AbstractTopDownSet",
          function(object, by=list(Mz=object$Mz, AgcTarget=object$AgcTarget)) {
    object@colData$MedianIonInjectionTimeMs <- tryCatch(
        ave(
            object@colData$IonInjectionTimeMs,
            .groupByLabels(by),
            FUN=function(x)median(x, na.rm=TRUE)
        ),
        warning=function(w) {
            stop("converted from warning: ", conditionMessage(w))
        }
    )
    object@colData$MedianIonInjectionTimeMs <- .colsToRle(
        object@colData[, "MedianIonInjectionTimeMs", drop=FALSE]
    )$MedianIonInjectionTimeMs
    msg <- "Recalculate median injection time"
    if (length(names(by))) {
        msg <- paste0(msg, " based on: ", paste0(names(by), collapse=", "))
    }
    .atdsLogMsg(object, msg, ".", addDim=FALSE)
})
