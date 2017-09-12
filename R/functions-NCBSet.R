#' Validate `NCBSet`
#'
#' @param object `NCBSet`
#' @return `TRUE` (if valid) else character with msg what was incorrect
#' @noRd
.validateNCBSet <- function(object) {
    msg <- character()

    if (!all(object@assay@x %in% as.double(0:3))) {
        msg <- c(msg, "'assay' contains invalid values.")
    }

    if (length(msg)) {
        msg
    } else {
        TRUE
    }
}
