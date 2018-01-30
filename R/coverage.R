#' Find best combination of columns for highest coverage
#'
#' @param x `dgCMatrix`
#' @param intensity `double` total intensity
#' @param n `integer`, max number of combinations/iterations
#' @param minN `integer` stop if there are less than `minN` non-zero
#' elements added by the next column.
#' @param maximise `character`, maximise fragment/bond coverage
#' @return `matrix`, first column: index, second column: number of fragments,
#' third column: number of bonds
#' @noRd
.bestNcbCoverageCombination <- function(x, intensity=NULL,
                                        n=ncol(x), minN=0L,
                                        maximise=c("fragments", "bonds")) {
    stopifnot(is(x, "dgCMatrix"))
    if (is.null(intensity)) {
        intensity <- rep(0L, ncol(x))
    }
    stopifnot(ncol(x) == length(intensity))
    maximise <- match.arg(maximise)
    m <- matrix(NA_real_, nrow=n, ncol=3L,
                dimnames=list(NULL, c("index", "fragments", "bonds")))

    ## make bonds binary but not logical because we want to use
    ## .removeNcbCombinations latter
    bonds <- x
    bonds@x[] <- 1L

    for (i in seq_len(n)) {
        if (maximise == "fragments") {
            hc <- .highestNcbFragmentCoverage(x, intensity=intensity)
            m[i, ] <- c(hc, .colCounts(bonds[, hc["index"], drop=FALSE]))
        } else {
            hc <- .highestNcbBondCoverage(bonds, intensity=intensity)
            m[i, c("index", "bonds", "fragments")] <-
                c(hc, .countFragments(x[, hc["index"], drop=FALSE]))
        }
        x <- .removeNcbCombinations(x, hc["index"])
        bonds <- .removeNcbCombinations(bonds, hc["index"])

        if (((!nnzero(x) && maximise == "fragments") ||
             (!nnzero(bonds) && maximise == "bonds")) &&
            hc[maximise] >= minN) {
            m <- m[seq_len(i), , drop=FALSE]
            break
        } else if (hc[maximise] < minN) {
            m <- m[seq_len(i - 1L), , drop=FALSE]
            break
        }
    }
    m
}

#' Highest NCB bond coverage
#'
#' Find column with highest coverage of bonds.
#' If multiple max. are found choose the one with the highest total intensity.
#'
#' @param x `dgCMatrix`
#' @param intensity `double` total intensity
#' @return `numeric`, first element: index of column with highest coverage,
#' second element: number bonds
#' @noRd
.highestNcbBondCoverage <- function(x, intensity=rep(0L, ncol(x))) {
    nb <- unname(.colCounts(x))

    ## we don't use `which.max` because if there are multiple matches we want to
    ## select the one with the highest intensity and not the first one.
    i <- which(nb == max(nb), useNames=FALSE)

    if (length(i) > 1L) {
        i <- i[which.max(intensity[i])]
    }
    c(index=i, bonds=nb[i])
}

#' Highest NCB fragment coverage
#'
#' Find column with highest coverage of NCB fragments, B fragments count twice.
#' If multiple max. are found choose the one with the highest total intensity.
#'
#' @param x `dgCMatrix`
#' @param intensity `double` total intensity
#' @param maximise `character`, maximise fragment/bond coverage
#' @return `numeric`, first element: index of column with highest coverage,
#' second element: number fragments
#' @noRd
.highestNcbFragmentCoverage <- function(x, intensity=rep(0L, ncol(x))) {
    nf <- unname(.countFragments(x))

    ## we don't use `which.max` because if there are multiple matches we want to
    ## select the one with the highest intensity and not the first one.
    i <- which(nf == max(nf), useNames=FALSE)

    if (length(i) > 1L) {
        i <- i[which.max(intensity[i])]
    }
    c(index=i, fragments=nf[i])
}

#' Remove NCB combinations
#'
#' @param x `dgCMatrix`
#' @param i `integer`, column with highest coverage
#' @return `dgCMatrix`, coverage reduced
#' @noRd
.removeNcbCombinations <- function(x, i) {
    stopifnot(is(x, "dgCMatrix"))
    stopifnot(!any(x@x > 3L))
    r <- .row(x[, i, drop=FALSE])
    hc <- x[r, i]
    r <- r - 1L
    ## if 3 (bidirectional) remove all
    x@x[x@i %in% r[hc == 3L]] <- 0L
    ## if 1 (N) keep 2 (C) and reduce 3 (bidirectional) to 2 (C)
    isN <- x@i %in% r[hc == 1L]
    isN <- isN & x@x != 2L
    x@x[isN] <- x@x[isN] - 1L
    ## if 2 (C) keep 1 (N) and reduce 3 (bidirectional) to 1 (N)
    isC <- x@i %in% r[hc == 2L]
    isC <- isC & x@x != 1L
    x@x[isC] <- x@x[isC] - 2L
    drop0(x)
}
