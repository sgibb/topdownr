#' mask matrix for matrix multiplication in .rowMeansGroup
#'
#' @param x character/numeric, group identifier
#' @return sparseMatrix with group masked
.createMaskMatrix <- function(x) {
  if (!is.numeric(x)) {
    x <- match(x, unique(x))
  }
  sparseMatrix(i=seq_along(x), j=x, x=1L)
}

#' rowMeans groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x Matrix
#' @param group group
#' @param ... further arguments passed to rowMeans
#' @return sparseMatrix
.rowMeansGroup <- function(x, group, ...) {
  stopifnot(is(x, "Matrix"))
  stopifnot(ncol(x) == length(group))
  mm <- .createMaskMatrix(group)
  m <- (x %*% mm)
  m@x <-m@x / ((x != 0L) %*% mm)@x
  m
}

