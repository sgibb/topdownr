#' Cumulative combination of conditions
#'
#' Used for visualization in the `fragmentationMap`
#'
#' @param x `dgcMatrix`
#' @return `dgcMatrix`
#' @noRd
.cumComb <- function(x) {
    stopifnot(is(x, "dgCMatrix") || is.matrix(x))
    m <- x
    x1 <- x == 1L
    x2 <- x == 2L
    ## use for loop instead of apply because the latter returns a vector instead
    ## of a matrix for a one-column input matrix (vs a matrix for a multi-column
    ## input, as.matrix/t just do different things for vectors/matrices).
    for (i in seq_len(nrow(x))) {
        m[i, ] <- cummax(x1[i, ]) + cummax(x2[i, ]) * 2L
    }
    m[x == 3L] <- 3L
    for (i in seq_len(nrow(m))) {
        m[i, ] <- cummax(m[i, ])
    }
    drop0(m)
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
