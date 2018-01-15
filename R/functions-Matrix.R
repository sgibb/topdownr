#' Column index
#'
#' @param x `dgCMatrix`
#' @return `numeric`, column index
#' @noRd
.col <- function(x) {
    stopifnot(is(x, "CsparseMatrix"))
    dp <- diff(x@p)
    rep(seq_along(dp), dp)
}

#' Count nonzero values in columns
#'
#' @param x `dgCMatrix`
#' @return `numeric`
#' @noRd
.colCounts <- function(x) {
    tabulate(.col(x), ncol(x))
}

#' colSums groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x `Matrix`
#' @param group `character`/`numeric`, group identifier
#' @return `sparseMatrix`
#' @noRd
.colSumsGroup <- function(x, group) {
    stopifnot(is(x, "Matrix"))
    stopifnot(nrow(x) == length(group))
    crossprod(.createMaskMatrix(group), x)
}

#' Count NCB fragments
#'
#' @param x `dgCMatrix`
#' @return `integer`, number of fragments per column
#' @noRd
.countFragments <- function(x) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(!any(x@x > 3L))
    i <- x@x > 1L
    x@x[i] <- x@x[i] - 1L
    Matrix::colSums(x)
}

#' mask matrix for matrix multiplication in .rowMeansGroup
#'
#' @param x `character`/`numeric`, group identifier
#' @return `sparseMatrix` with group masked
#' @noRd
.createMaskMatrix <- function(x) {
    if (!is.numeric(x)) {
        x <- match(x, unique(x))
    }
    sparseMatrix(i=seq_along(x), j=x, x=1L)
}

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

#' convert sparseMatrix to data.frame (for ggplot2)
#'
#' @param x `dgCMatrix`
#' @return `data.frame`
#' @noRd
.dgcMatrix2data.frame <- function(x) {
    data.frame(row=.row(x), col=.col(x), x=x@x, stringsAsFactors=FALSE)
}

#' drop0 but row-wise with length(tol) == nrow(x)
#'
#' @param x `dgCMatrix`
#' @param tol `sparseVector`/`double`, length==nrow(x), less than
#' @return `dgCMatrix`
#' @noRd
.drop0rowLe <- function(x, tol) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(length(tol) == nrow(x))
    x@x[x@x <= as.vector(tol)[.row(x)]] <- 0L
    drop0(x, tol=0L, is.Csparse=TRUE)
}

.drop0rowLt <- function(x, tol) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(length(tol) == nrow(x))
    x@x[x@x < as.vector(tol)[.row(x)]] <- 0L
    drop0(x, tol=0L, is.Csparse=TRUE)
}

#' drop0 but just row-wise groups with fewer than minN entries per group
#'
#' @param x `dgCMatrix`
#' @param group `character`/`numeric`, group identifier
#' @param minN `integer`, minimal number of group members > 0
#' @return `dgCMatrix`
#' @noRd
.drop0rowReplicates <- function(x, group, minN) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(length(group) == ncol(x))
    stopifnot(is.integer(minN) && length(minN) == 1L)
    m <- .createMaskMatrix(group)
    l <- x
    l@x[] <- 1L
    l <- tcrossprod(l %*% m, # identical to .rowSumsGroup
                    m)
    l@x[l@x < minN] <- 0L
    l <- drop0(l, tol=0L, is.Csparse=TRUE)
    l@x[] <- 1L
    x * l
}

#' drop0 but for `NA`
#'
#' @param x `dgCMatrix`
#' @return `dgCMatrix`
#' @noRd
.dropNA <- function(x) {
    stopifnot(is(x, "dgCMatrix"))
    x@x[is.na(x@x)] <- 0L
    drop0(x, tol=0L, is.Csparse=TRUE)
}

#' normalise (col-wise scale)
#'
#' @param x `dgCMatrix`
#' @param scale `double`, scale variable
#' @return `dgCMatrix`
#' @noRd
.normaliseCols <- function(x, scale=.rowMax(t(x))) {
    t(.normaliseRows(t(x), scale=scale))
}

#' normalise (row-wise scale)
#'
#' @param x `dgCMatrix`
#' @param scale `double`, scale variable
#' @return `dgCMatrix`
#' @noRd
.normaliseRows <- function(x, scale=.rowMax(x)) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(
        (is.numeric(scale) || is(scale, "sparseVector")) &&
        (length(scale) == 1L || length(scale) == nrow(x))
    )
    x@x <- x@x / scale[.row(x)]
    x
}

#' Row index
#'
#' @param x `dgCMatrix`
#' @return `numeric`, row index
#' @noRd
.row <- function(x) {
    stopifnot(is(x, "CsparseMatrix"))
    x@i + 1L
}

#' Count nonzero values in rows
#'
#' @param x `dgCMatrix`
#' @return `numeric`
#' @noRd
.rowCounts <- function(x) {
    tabulate(.row(x), nbins=nrow(x))
}

#' rowCvs groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x `Matrix`
#' @param group `integer`/`character` group identifier
#' @param na.rm `logical`, should `NA`s removed?
#' @return `sparseMatrix`
#' @noRd
.rowCvsGroup <- function(x, group, na.rm=TRUE) {
    y <- .rowSdsGroup(x, group, na.rm)
    i <- Matrix::which(y != 0L)
    m <- .rowMeansGroup(x, group, na.rm)
    y@x <- y@x / m[i]
    y
}

#' rowMax row-wise maximum, similar to apply(x, 1, max) but faster on
#' sparseMatrix
#'
#' @param x `dgCMatrix`
#' @param na.rm `logical`, should `NA`s removed?
#' @return `sparseVector`
#' @noRd
.rowMax <- function(x, na.rm=TRUE) {
    stopifnot(is(x, "dgCMatrix"))
    sparseVector(
        .vapply1d(split(x@x, .row(x)), max, na.rm=na.rm),
        i=sort.int(unique(.row(x))), length=nrow(x)
    )
}

#' rowMeans groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x `Matrix`
#' @param group `integer`/`character` group identifier
#' @param na.rm `logical`, should `NA`s removed?
#' @return `sparseMatrix`
#' @noRd
.rowMeansGroup <- function(x, group, na.rm=TRUE) {
    stopifnot(is(x, "Matrix"))
    stopifnot(ncol(x) == length(group))
    if (na.rm) {
        x <- .dropNA(x)
    }
    mm <- .createMaskMatrix(group)
    m <- x %*% mm
    m@x <- m@x / ((x != 0L) %*% mm)@x
    m
}

#' rowSds groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x `Matrix`
#' @param group `integer`/`character` group identifier
#' @param na.rm `logical`, should `NA`s removed?
#' @return `sparseMatrix`
#' @noRd
.rowSdsGroup <- function(x, group, na.rm=TRUE) {
    stopifnot(is(x, "Matrix"))
    stopifnot(ncol(x) == length(group))

    if (na.rm) {
        nna <- .dropNA(x)
    } else {
        nna <- x
    }
    nna@x[] <- 1L
    nna <- .rowSumsGroup(nna, group=group, na.rm=na.rm)
    nna@x[nna@x < 2L] <- NA_real_ # return NA if n < 2 (similar to sd)
    var <- .rowMeansGroup(x * x, group=group, na.rm=na.rm) -
        .rowMeansGroup(x, group=group, na.rm=na.rm)^2L
    nna1 <- nna
    nna1@x <- nna1@x - 1L
    var@x <- (var * nna)@x / nna1@x
    # we want a sparse matrix; NA == missing, we dropNA
    .dropNA(sqrt(var))
}

#' rowSums groupwise, similar to rowsum but for sparceMatrices
#'
#' @param x `Matrix`
#' @param group `integer`/`character` group identifier
#' @param na.rm `logical`, should `NA`s removed?
#' @return `sparseMatrix`
#' @noRd
.rowSumsGroup <- function(x, group, na.rm=TRUE) {
    stopifnot(is(x, "Matrix"))
    stopifnot(ncol(x) == length(group))
    if (na.rm) {
        x <- .dropNA(x)
    }
    x %*% .createMaskMatrix(group)
}

#' summary
#'
#' @param x `Matrix`
#' @param what `character`, "rows"/"columns"
#' @param na.rm `logical`, should NA removed?
#' @return `data.frame`
#' @noRd
.summary <- function(x, what=c("rows", "columns"), na.rm=TRUE) {

    if (match.arg(what) == "columns") {
        x <- t(x)
    }

    i <- sort.int(unique(.row(x)))
    qq <- matrix(0, nrow=5L, ncol=nrow(x))
    qq[, i] <- vapply(
        X=split(x@x, .row(x)), FUN=quantile, FUN.VALUE=double(5L),
        probs=seq(0, 1, 0.25), na.rm=na.rm, USE.NAMES=FALSE
    )

    data.frame(
        Fragments=.rowCounts(x),
        Total=Matrix::rowSums(x, na.rm=na.rm),
        Min=qq[1L, ],
        Q1=qq[2L, ],
        Median=qq[3L, ],
        Mean=Matrix::rowMeans(x, na.rm=na.rm),
        Q3=qq[4L, ],
        Max=qq[5L, ]
    )
}
