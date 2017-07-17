#' aggregate data.frame (combine numeric columns by \code{fun}, and take the
#' first row for non-numeric columns
#'
#' @param x data.frame
#' @param f character, grouping value
#' @param ignoreNumCols character, column names that won't be aggregated by
#' \code{fun}
#' @param fun function, aggregation function
#' @return data.frame
#' @noRd
.aggregateDataFrame <- function(x, f, ignoreNumCols=character(), fun=mean, na.rm=TRUE) {
  fun <- match.fun(fun)

  cn <- colnames(x)

  isNumCol <- .vapply1l(x, function(column) {
     is.numeric(column) || (is(column, "Rle") && is.numeric(runValue(column)))
  }) & !cn %in% ignoreNumCols

  nonNum <- x[!duplicated(f), !isNumCol, drop=FALSE]
  rn <- rownames(nonNum)
  num <- aggregate(x[, isNumCol, drop=FALSE], by=list(f), FUN=fun, na.rm=na.rm,
                   drop=FALSE)
  x <- .colsToRle(cbind(nonNum, num)[, cn])
  rownames(x) <- rn
  x
}

#' Convert DataFrame columns to Rle
#'
#' @param x DataFrame
#' @return DataFrame
#' @noRd
.colsToRle <- function(x) {
  toConvert <- .vapply1l(x, function(xx)length(unique(xx)) < nrow(x) / 4L)
  x[toConvert] <- lapply(x[toConvert], Rle)
  x
}

#' droplevels for Rle/factor columns
#'
#' @param x DataFrame
#' @return DataFrame
#' @noRd
.droplevels <- function(x) {
  isFactorColumn <- .vapply1l(x, function(column) {
     is.factor(column) || (is(column, "Rle") && is.factor(runValue(column)))
  })
  x[isFactorColumn] <- droplevels(x[isFactorColumn])
  x
}
