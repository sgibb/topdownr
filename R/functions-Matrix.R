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

#' drop0 but row-wise with length(tol) == nrow(x)
#'
#' @param x dgCMatrix
#' @param tol sparseVector/double, length==nrow(x), less than
#' @return dgCMatrix
#' @noRd
.drop0rowLe <- function(x, tol) {
  stopifnot(is(x, "dgCMatrix"))
  stopifnot(length(tol) == nrow(x))
  x@x[x@x <= as.vector(tol)[x@i + 1L]] <- 0L
  drop0(x, tol=0L, is.Csparse=TRUE)
}

.drop0rowLt <- function(x, tol) {
  stopifnot(is(x, "dgCMatrix"))
  stopifnot(length(tol) == nrow(x))
  x@x[x@x < as.vector(tol)[x@i + 1L]] <- 0L
  drop0(x, tol=0L, is.Csparse=TRUE)
}

#' calculate rect coordinates for a Matrix
#'
#' @param x Matrix
#' @param width double
#' @param height double
#' @return matrix with coords (xleft, ybottom, xright, ytop) (row index is x,
#' col index is y), and color
#' @noRd
.m2rect <- function(x, width=1L, height=1L) {
  stopifnot(is(x, "Matrix"))
  dp <- diff(x@p)
  y <- rep(seq_along(dp), dp) - 1L
  cbind(xleft=x@i, ybottom=y, xright=x@i + width, ytop=y + height, col=x@x)
}

#' normalise (row-wise scale) to 0:1
#'
#' @param x `dgCMatrix`
#' @return `dgCMatrix`
#' @noRd
.normaliseRows <- function(x) {
  stopifnot(is(x, "dgCMatrix"))
  scale <- .rowMax(x)
  x@x <- x@x / scale[x@i + 1L]
  x
}

#' rowMax row-wise maximum, similar to apply(x, 1, max) but faster on
#' sparseMatrix
#'
#' @param x dgCMatrix
#' @return sparseVector
#' @noRd
.rowMax <- function(x) {
  stopifnot(is(x, "dgCMatrix"))
  sparseVector(.vapply1d(split(x@x, x@i), max),
               i=sort.int(unique(x@i)) + 1L, length=nrow(x))
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

