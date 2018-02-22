#' Add log message.
#'
#' @param object `AbstractTopDownSet`
#' @param addDim `logical`, add dimension information
#' @return `AbstractTopDownSet`
#' @noRd
.atdsLogMsg <- function(object, ..., addDim=TRUE) {
    if (!isTRUE(inherits(object, "AbstractTopDownSet"))) {
        stop("'object' has to be an 'AbstractTopDownSet' object.")
    }
    if (addDim) {
        object@processing <- c(
            object@processing,
            .logmsg(..., "; ", .logdim(object), "."))
    } else {
        object@processing <- c(object@processing, .logmsg(...))
    }
    object
}

#' Get fragment types
#'
#' @param object `AbstractTopDownSet`
#' @return `character`
#' @noRd
fragmentType <- function(object) {
    .inheritsAbstractTopDownSet(object)
    elementMetadata(object@rowViews)$type
}

#' Test for AbstractTopDownSetclass
#'
#' @param object object to test
#' @return `TRUE` if object is an AbstractTopDownSet otherwise fails with an
#' error
#' @noRd
.inheritsAbstractTopDownSet <- function(object) {
    if (!isTRUE(inherits(object, "AbstractTopDownSet"))) {
        stop("'object' doesn't inherit 'AbstractTopDownSet'.")
    }
    TRUE
}

#' String with number of fragments and dimension of the assay.
#'
#' @param object `TopDownSet`
#' @return `character`
#' @noRd
.logdim <- function(object) {
    .inheritsAbstractTopDownSet(object)
    sprintf("%d fragments [%d;%d]",
            nnzero(object@assay), nrow(object), ncol(object))
}

#' Validate `AbstractTopDownSet`
#'
#' @param object `AbstractTopDownSet`
#' @return `TRUE` (if valid) else character with msg what was incorrect
#' @noRd
.validateAbstractTopDownSet <- function(object) {
    msg <- character()

    if (nrow(object@assay) != length(object@rowViews)) {
        msg <- c(
            msg,
            "Mismatch between fragment data in 'rowViews' and 'assay'."
        )
    }

    if (any(rownames(object@assay) != names(object@rowViews))) {
        msg <- c(
            msg,
            "Mismatch between fragment names in 'rowViews' and 'assay'."
        )
    }

    if (ncol(object@assay) != nrow(object@colData)) {
        msg <- c(
            msg,
            "Mismatch between condition data in 'colData' and 'assay'."
        )
    }

    if (any(colnames(object@assay) != rownames(object@colData))) {
        msg <- c(
            msg,
            "Mismatch between condition names in 'colData' and 'assay'."
        )
    }

    if (length(msg)) {
        msg
    } else {
        TRUE
    }
}
