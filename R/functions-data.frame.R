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
  rn <- rownames(x)

  isNumCol <- .vapply1l(x, is.numeric) & !cn %in% ignoreNumCols

  rn <- .vapply1c(split(rn, f), "[", 1L)
  nonNum <- aggregate(x[, !isNumCol, drop=FALSE], by=f, FUN="[[", 1L,
                      drop=FALSE, simplify=TRUE)
  num <- aggregate(x[, isNumCol, drop=FALSE], by=f, FUN=fun, na.rm=na.rm,
                   drop=FALSE, simplify=TRUE)
  cbind(nonNum, num, stringsAsFactors=FALSE, row.names=rn)[, cn]
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
