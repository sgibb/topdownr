#' rowMeans groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x Matrix
#' @param group group
#' @param ... further arguments passed to rowMeans
#' @return sparseMatrix
.rowMeansGroup <- function(x, group, ...) {
  stopifnot(ncol(x) == length(group))
  j <- seq_len(ncol(x))
  l <- lapply(split(j, group), function(jj) {
    as(Matrix::rowMeans(x[, jj, drop=FALSE], na.rm=TRUE, sparseResult=TRUE),
       "sparseMatrix")
  })
  do.call(cbind, l)
}
