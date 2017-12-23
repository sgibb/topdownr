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
    x@files <- x@files[.subsetFiles(x@files, unique(x@colData$File)) | isFasta]

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
