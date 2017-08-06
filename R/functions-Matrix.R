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

#' drop0 but row-wise with length(tol) == nrow(x)
#'
#' @param x `dgCMatrix`
#' @param tol `sparseVector`/`double`, length==nrow(x), less than
#' @return `dgCMatrix`
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

#' calculate rect coordinates for a Matrix
#'
#' @param x `Matrix`
#' @param width `double`
#' @param height `double`
#' @return `matrix` with coords (xleft, ybottom, xright, ytop) (row index is x,
#' col index is y), and color
#' @noRd
.m2rect <- function(x, width=1L, height=1L) {
    stopifnot(is(x, "Matrix"))
    dp <- diff(x@p)
    y <- rep(seq_along(dp), dp) - 1L
    cbind(xleft=x@i, ybottom=y, xright=x@i + width, ytop=y + height, col=x@x)
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
    stopifnot((is.numeric(scale) || is(scale, "sparseVector")) &&
              (length(scale) == 1L || length(scale) == nrow(x)))
    x@x <- x@x / scale[x@i + 1L]
    x
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
    sparseVector(.vapply1d(split(x@x, x@i), max, na.rm=na.rm),
                 i=sort.int(unique(x@i)) + 1L, length=nrow(x))
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
    m <- (x %*% mm)
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
