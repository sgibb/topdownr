#' Cumulative combination of conditions
#'
#' Used for visualization in the `fragmentationMap`
#'
#' @param x `dgcMatrix`
#' @return `dgcMatrix`
#' @noRd
.cumComb <- function(x) {
    stopifnot(is(x, "dgCMatrix") || is.matrix(x))
    x <- as.matrix(x)
    m <- t(apply(x == 1L, 1L, cummax) + apply(x == 2L, 1L, cummax) * 2L)
    m[x == 3L] <- 3L
    drop0(t(apply(m, 1L, cummax)))
}

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
