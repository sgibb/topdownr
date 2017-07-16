#' aggregate data.frame (combine numeric columns by \code{fun}, and take the
#' first row for non-numeric columns
#' @param x data.frame
#' @param f character, grouping value
#' @param ignoreCols character, column names that won't be part of the returned
#' data.frame
#' @return data.frame
#' @noRd
.aggregateDataFrame <- function(x, f, ignoreCols=character(), na.rm=TRUE) {

  x <- x[, !colnames(x) %in% ignoreCols]

  cn <- colnames(x)
  isNumCol <- .vapply1l(x, is.numeric)
  nonDuplicatedF <- !duplicated(f)
  rn <- rownames(x)[nonDuplicatedF]
  num <- MSnbase:::rowmean(x[, isNumCol, drop=FALSE], group=f, reorder=FALSE,
                           na.rm=na.rm)
  nonNum <- x[nonDuplicatedF, !isNumCol, drop=FALSE]
  x <- cbind(num, nonNum, stringsAsFactors=FALSE, row.names=rn)[, cn]
}
