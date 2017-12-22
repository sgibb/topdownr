#' aggregate data.frame (combine numeric columns by fun, and take the
#' first row for non-numeric columns
#'
#' @param x `data.frame`
#' @param f `character`, grouping value
#' @param ignoreNumCols `character`, column names that won't be aggregated by
#' fun
#' @param fun `function`, aggregation function
#' @return `data.frame`
#' @noRd
.aggregateDataFrame <- function(x, f, ignoreNumCols=character(),
                                fun=mean, na.rm=TRUE) {
    fun <- match.fun(fun)

    cn <- colnames(x)

    isNumCol <- .isNumCol(x) & !cn %in% ignoreNumCols

    nonNum <- x[!duplicated(f), !isNumCol, drop=FALSE]
    rn <- rownames(nonNum)
    num <- aggregate(
        x[, isNumCol, drop=FALSE], by=list(f), FUN=fun, na.rm=na.rm, drop=FALSE
    )
    ## resort (aggregate turns "by" into a factor (locale depended sorting))
    num <- num[match(unique(f), num[["Group.1"]]),, drop=FALSE]
    x <- .colsToRle(cbind(nonNum, num)[, cn])
    rownames(x) <- rn
    x
}

#' Convert DataFrame columns to logical
#'
#' @param x `DataFrame`
#' @return `DataFrame`
#' @noRd
.colsToLogical <- function(x) {
    toConvert <- .isCharacterCol(x)
    x[toConvert] <- lapply(x[toConvert], .characterToLogical)
    x
}

#' Convert DataFrame columns to Rle
#'
#' @param x `DataFrame`
#' @return `DataFrame`
#' @noRd
.colsToRle <- function(x) {
    toConvert <- .vapply1l(x, function(xx)length(unique(xx)) < nrow(x) / 4L)
    x[toConvert] <- lapply(x[toConvert], Rle)
    x
}

#' droplevels for Rle/factor columns
#'
#' @param x `DataFrame`
#' @return `DataFrame`
#' @noRd
.droplevels <- function(x) {
    isFactorColumn <- .vapply1l(x, function(column) {
        is.factor(column) || (is(column, "Rle") && is.factor(runValue(column)))
    })
    x[isFactorColumn] <- droplevels(x[isFactorColumn])
    x
}

#' Drop non informative columns (all rows are identical)
#'
#' @param x `data.frame`/`DataFrame`
#' @param keep `character` column names that should never be dropped.
#' @return x, without columns that are identical
#' @noRd
.dropNonInformativeColumns <- function(x, keep="Mz") {
    keep <- !.vapply1l(x, .allIdentical) | colnames(x) %in% keep
    x[, keep, drop=FALSE]
}

#' Test for character columns
#'
#' @param x `data.frame`
#' @return `logical`
#' @noRd
.isCharacterCol <- function(x) {
    .vapply1l(x, function(column) {
        is.character(column) ||
            (is(column, "Rle") && is.character(runValue(column)))
    })
}

#' Test for numeric columns
#'
#' @param x `data.frame`
#' @return `logical`
#' @noRd
.isNumCol <- function(x) {
    .vapply1l(x, function(column) {
        is.numeric(column) ||
            (is(column, "Rle") && is.numeric(runValue(column)))
    })
}

#' Order data.frame by multiple columns given as character vector
#'
#' @param x `data.frame`
#' @param cols `character`, column names
#' @return `integer`
#' @noRd
.orderByColumns <- function(x, cols) {
    stopifnot(is.data.frame(x) || inherits(x, "DataFrame"))
    stopifnot(all(cols %in% colnames(x)))
    do.call(order, x[cols])
}

#' Combine data.frames rowwise, similar to base::rbind but creates missing
#' columns
#'
#' @param ... `data.frame`s
#' @noRd
.rbind <- function(...) {
    l <- list(...)

    if (length(l) == 1L) {
        l <- l[[1L]]
    }

    stopifnot(
        all(.vapply1l(l, function(ll) {
            is.data.frame(ll) || inherits(ll, "DataFrame")
        }))
    )

    nms <- lapply(l, names)
    allcn <- unique(unlist(nms))

    for (i in seq(along=l)) {
        cnMissing <- setdiff(allcn, nms[[i]])
        if (length(cnMissing)) {
            l[[i]][, cnMissing] <- NA
        }
    }
    do.call(rbind, l)
}
