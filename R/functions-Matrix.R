#' mask matrix for matrix multiplication in .rowMeansGroup
#'
#' @param x character/numeric, group identifier
#' @return sparseMatrix with group masked
#' @noRd
.createMaskMatrix <- function(x) {
  if (!is.numeric(x)) {
    x <- match(x, unique(x))
  }
  sparseMatrix(i=seq_along(x), j=x, x=1L)
}

#' colSums groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x Matrix
#' @param group group
#' @param ... further arguments passed to rowMeans
#' @return sparseMatrix
#' @noRd
.colSumsGroup <- function(x, group) {
  stopifnot(is(x, "Matrix"))
  stopifnot(nrow(x) == length(group))
  mm <- .createMaskMatrix(group)
  t(mm) %*% x
}

#' calculate rect coordinates for a Matrix
#'
#' @param x Matrix
#' @param width double
#' @param height double
#' @return matrix with coords (xleft, ybottom, xright, ytop) (row index is x,
#' col index is y), and color
.m2rect <- function(x, width=1L, height=1L) {
  stopifnot(is(x, "Matrix"))
  dp <- diff(x@p)
  y <- rep(seq_along(dp), dp) - 1L
  w2 <- width / 2L
  h2 <- height / 2L
  cbind(xleft=x@i - w2, ybottom=y - h2, xright=x@i + w2, ytop=y + h2, col=x@x)
}

#' rowMeans groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x Matrix
#' @param group group
#' @return sparseMatrix
#' @noRd
.rowMeansGroup <- function(x, group) {
  stopifnot(is(x, "Matrix"))
  stopifnot(ncol(x) == length(group))
  mm <- .createMaskMatrix(group)
  m <- (x %*% mm)
  m@x <-m@x / ((x != 0L) %*% mm)@x
  m
}

