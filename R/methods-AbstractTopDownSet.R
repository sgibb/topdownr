#' @param x `AbstractTopDownSet`
#' @param i `numeric`, `logical` or `character`, subsetting on fragment data,
#' names (`c("a1", "b1", "c1", "c2", "c3")`) and types (`c("c", "x")`) are
#' supported.
#' @param j `numeric` or `logical`, subsetting based on condition data.
#' @param \ldots currently ignored.
#' @param drop `logical`, currently ignored.
#' @noRd
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

    x <- .atdsLogMsg(x, "Subsetted ", ld0, " to ", ld1, ".", addDim=FALSE)

    if (validObject(x)) {
        x
    }
})

#' @param x `AbstractTopDownSet`
#' @param i `numeric`, `logical` or `character`, subsetting based on condition
#' data.
#' @param j currently ignored.
#' @param \ldots currently ignored.
#' @noRd
setMethod("[[", c("AbstractTopDownSet", "ANY", "missing"),
          function(x, i, j, ...) {
    colData(x)[[i, ...]]
})

#' @param x [AbstractTopDownSet-class]
#' @noRd
setReplaceMethod("[[", c("AbstractTopDownSet", "ANY", "missing"),
                 function(x, i, j, ..., value) {
    colData(x)[[i, ...]] <- value
    if (validObject(x)) {
      x
    }
})

#' @param x `AbstractTopDownSet`
#' @noRd
#' @export
.DollarNames.AbstractTopDownSet <- function(x, pattern="") {
    grep(pattern, names(colData(x)), value=TRUE)
}

#' @param x `AbstractTopDownSet`
#' @noRd
setMethod("$", "AbstractTopDownSet", function(x, name) {
    colData(x)[[name]]
})

#' @param x `AbstractTopDownSet`
#' @noRd
setReplaceMethod("$", "AbstractTopDownSet", function(x, name, value) {
    colData(x)[[name]] <- value
    x
})

#' @param object `AbstractTopDownSet`
#' @return `Matrix`
#' @export
#' @noRd
setMethod("assayData", "AbstractTopDownSet", function(object) {
    object@assay
})

#' @param object `AbstractTopDownSet`
#' @return `DataFrame`
#' @export
#' @noRd
setMethod("colData", "AbstractTopDownSet", function(object) {
    object@colData
})

#' @param x `AbstractTopDownSet`
#' @return `AbstractTopDownSet`
#' @export
#' @noRd
setReplaceMethod("colData", "AbstractTopDownSet", function(object, ..., value) {
    object@colData <- value
    if (validObject(object)) {
      object
    }
})

#' @param object `AbstractTopDownSet`
#' @return `DataFrame`
#' @export
#' @noRd
setMethod("conditionData", "AbstractTopDownSet", function(object, ...) {
    colData(object)
})

#' @param object `AbstractTopDownSet`
#' @return `AbstractTopDownSet`
#' @export
#' @noRd
setReplaceMethod("conditionData", "AbstractTopDownSet",
                 function(object, ..., value) {
    colData(object) <- value
    if (validObject(object)) {
      object
    }
})

#' @param object `AbstractTopDownSet`
#' @return `numeric`
#' @noRd
setMethod("dim", "AbstractTopDownSet", function(x) {
    dim(x@assay)
})

#' @param object `AbstractTopDownSet`
#' @return `list`
#' @noRd
setMethod("dimnames", "AbstractTopDownSet", function(x) {
    list(names(x@rowViews), row.names(x@colData))
})

#' @param object `AbstractTopDownSet`
#' @return `XStringViews`
#' @export
#' @noRd
setMethod("fragmentData", "AbstractTopDownSet", function(object, ...) {
    rowViews(object)
})

#' @param object `AbstractTopDownSet`
#' @return `AbstractTopDownSet`
#' @export
#' @noRd
setMethod("removeEmptyConditions", "AbstractTopDownSet",
          function(object) {
    i <- Matrix::colSums(object@assay) != 0L
    object <- object[, i]
    .atdsLogMsg(object, sum(!i), " empty conditions removed")
})

#' @param object `AbstractTopDownSet`
#' @return `XStringViews`
#' @export
#' @noRd
setMethod("rowViews", "AbstractTopDownSet", function(object, ...) {
    object@rowViews
})

#' @param object `AbstractTopDownSet`
#' @noRd
setMethod("show", "AbstractTopDownSet", function(object) {
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

    invisible(NULL)
})
