context("Matrix")

m <- sparseMatrix(i=rep(1:4, each=5), j=rep(1:10, 2), x=1:20)

test_that(".cbind", {
    expect_error(.cbind(m, matrix(1:10, ncol=2)))
    expect_error(.cbind(m, m))
    m1 <- sparseMatrix(i=1:4, j=1:4, x=1:4,
                       dimnames=list(LETTERS[1:4], NULL))
    expect_equal(.cbind(m1, m1), cbind(m1, m1))
    m2 <- sparseMatrix(i=1:4, j=1:4, x=c(1:2, 5:6),
                       dimnames=list(LETTERS[c(1:2, 5:6)], NULL))
    mr <- sparseMatrix(i=c(1:4, 1:2, 5:6), j=1:8, x=c(1:4, 1:2, 5:6),
                       dimnames=list(LETTERS[1:6], NULL))
    expect_equal(.cbind(m1, m2), mr)
})

test_that(".col", {
    expect_error(.col(matrix(1:10, ncol=2)))
    expect_equal(.col(m), rep(1:10, each=2))
})

test_that(".colCounts", {
    expect_error(.colCounts(matrix(1:10, ncol=2)))
    expect_equal(.colCounts(m), rep(2, 10))
})

test_that(".colSumsGroup", {
    r <- sparseMatrix(i=rep(1:2, each=5),
                      j=1:10,
                      x=seq(12, 30, by=2))
    r2 <- sparseMatrix(i=rep(1:2, 2),
                       j=1:4,
                       x=c(15, 40, 65, 90))
    expect_error(.colSumsGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(.colSumsGroup(m, group=1:2))
    expect_equal(.colSumsGroup(m, group=rep(1:2, 2)), r)
    expect_equal(.colSumsGroup(t(m), group=rep(1:2, each=5)), r2)
})

test_that(".countFragments", {
    expect_error(.countFragments(matrix(1:10, ncol=2)))
    expect_error(.countFragments(m))
    expect_equal(.countFragments(drop0(m %% 4)), rep(c(3, 1), 5))
})

test_that(".createMaskMatrix", {
    r1 <- sparseMatrix(i=1:10, j=rep(1:5, 2), x=1)
    r2 <- sparseMatrix(i=1:10, j=rep(1:2, each=5), x=1)
    expect_equal(.createMaskMatrix(rep(1:5, 2)), r1)
    expect_equal(.createMaskMatrix(rep(letters[1:5], 2)), r1)
    expect_equal(.createMaskMatrix(rep(1:2, each=5)), r2)
})

test_that(".cumComb", {
    a <- sparseMatrix(i=rep(1:3, 3:1),
                      j=c(1:2, 4, 3:4, 1),
                      x=c(1, 2, 2, 3, 1, 2))
    r <- sparseMatrix(i=rep(1:3, c(4, 2, 4)),
                      j=c(1:4, 3:4, 1:4),
                      x=rep(c(1, 3, 2), c(1, 5, 4)))
    a2 <- sparseMatrix(i=c(1:3),
                       j=rep(1, 3),
                       x=c(1, 2, 3))
    expect_error(.cumComb(1:10), "is not of class 'dgCMatrix'")
    expect_equal(.cumComb(a), r)
    expect_equal(.cumComb(a2), a2)
})

test_that(".dgCMatrix2data.frame", {
    d <- data.frame(row=c(rep(c(1, 3), 5), rep(c(2, 4), 5)),
                    col=rep(1:10, each=2), x=rep(1:10, each=2) + c(0, 10),
                    stringsAsFactors=FALSE)
    expect_error(.dgCMatrix2data.frame(matrix(1:10, nrow=2)))
    expect_equal(.dgCMatrix2data.frame(m), d)
})

test_that(".drop0rowLe/Lt", {
    r1 <- sparseMatrix(i=rep(1:4, each=3), j=rep(c(3:5, 8:10), 2),
                       x=c(3:5, 8:10, 13:15, 18:20))
    r2 <- sparseMatrix(i=1:10, j=rep(3:4, each=5), x=11:20)
    expect_error(.drop0rowLe(matrix(1:10, ncol=2), 1))
    expect_error(.drop0rowLt(matrix(1:10, ncol=2), 1))
    expect_error(.drop0rowLe(m, 1))
    expect_error(.drop0rowLt(m, 1))
    expect_equal(.drop0rowLe(m, c(2, 7, 12, 17)), r1)
    expect_equal(.drop0rowLt(m, c(3, 8, 13, 18)), r1)
    expect_equal(.drop0rowLe(t(m), 1:10), r2)
    expect_equal(.drop0rowLt(t(m), 2:11), r2)
})

test_that(".drop0rowReplicates", {
    group <- rep(1:5, each=2)
    expect_error(.drop0rowReplicates(matrix(1:10, ncol=2), 1))
    expect_error(.drop0rowReplicates(m, 1))
    expect_error(.drop0rowReplicates(m, group, 1.5))
    expect_error(.drop0rowReplicates(m, group, 1L:2L))
    expect_equal(.drop0rowReplicates(m, group, minN=1L), m)
    expect_equal(.drop0rowReplicates(m, group, minN=2L),
                 sparseMatrix(i=rep(1:4, each=4),
                              j=rep(c(1:4, 7:10), 2),
                              x=c(1:4, 7:10, 11:14, 17:20)))
    expect_equal(.drop0rowReplicates(m, group, minN=3L),
                 Matrix(0L, nrow=4, ncol=10, sparse=TRUE))

    expect_equal(.drop0rowReplicates(
                    sparseMatrix(i=c(1:2, 1:3, 1, 3, 2),
                                 j=rep(1:5, c(2, 3, 1, 1, 1)),
                                 x=2:9),
                    rep(1:2, 3:2), 2L),
                 sparseMatrix(i=c(1:2, 1:2, 1),
                              j=rep(1:3, c(2, 2, 1)),
                              x=c(2:5, 7), dims=c(3, 5)))
})

test_that(".dropNA", {
    r <- sparseMatrix(i=rep(1:4, each=3), j=rep(c(3:5, 8:10), 2),
                      x=c(3:5, 8:10, 13:15, 18:20))
    n <- m
    n[cbind(rep(1:4, each=2), rep(c(1:2, 6:7), 2))] <- NA
    expect_error(.dropNA(matrix(1:10, ncol=2)))
    expect_equal(.dropNA(m), m)
    expect_equal(.dropNA(n), r)
})

test_that(".normaliseCols", {
    expect_error(.normaliseCols(matrix(1:10, nrow=2)))
    expect_error(.normaliseCols(m, "A"))
    expect_error(.normaliseCols(m, 1:2))
    expect_equal(.normaliseCols(m),
                 as((scale(m, center=FALSE, scale=11:20)),
                    "dgCMatrix"))
    expect_equal(.normaliseCols(m, 11:20),
                 as((scale(m, center=FALSE, scale=11:20)),
                    "dgCMatrix"))
})

test_that(".normaliseRows", {
    expect_error(.normaliseRows(matrix(1:10, nrow=2)))
    expect_error(.normaliseRows(m, "A"))
    expect_error(.normaliseRows(m, 1:2))
    expect_equal(.normaliseRows(m),
                 as((t(scale(t(m), center=FALSE, scale=c(5, 10, 15, 20)))),
                    "dgCMatrix"))
    expect_equal(.normaliseRows(m, 1:4),
                 as((t(scale(t(m), center=FALSE, scale=1:4))),
                    "dgCMatrix"))
})

test_that(".row", {
    expect_error(.row(matrix(1:10, ncol=2)))
    expect_equal(.row(m), c(rep(c(1, 3), 5), rep(c(2, 4), 5)))
})

test_that(".rowCounts", {
    expect_error(.rowCounts(matrix(1:10, ncol=2)))
    expect_equal(.rowCounts(m), rep(5, 4))
})

test_that(".rowCvsGroup", {
    .rowcvs <- function(x, group, na.rm=TRUE) {
        l <- lapply(split(1L:ncol(x), group), function(i) {
            apply(x[, i, drop=FALSE], 1, function(xx) {
                      sd(xx, na.rm=na.rm) / mean(xx, na.rm=na.rm)
            })
        })
        do.call(cbind, l)
    }
    bm <- as.matrix(m)
    bm[bm == 0L] <- NA
    group <- rep(1:5, each=2)
    bmr <- .rowcvs(bm, group)
    dimnames(bmr) <- NULL
    mr <- as.matrix(.rowCvsGroup(m, group))
    mr[mr == 0L] <- NA
    dimnames(mr) <- NULL
    expect_equal(bmr, mr)

    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    bn <- as.matrix(n)
    bn[bn == 0L] <- NA
    group <- rep(1:2, each=5)
    bnr <- .rowcvs(bn, group, na.rm=FALSE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(.rowCvsGroup(n, group, na.rm=FALSE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)

    bnr <- .rowcvs(bn, group, na.rm=TRUE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(.rowCvsGroup(n, group, na.rm=TRUE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)
})

test_that(".rowMax", {
    expect_error(.rowMax(matrix(1:10, ncol=2)))
    expect_equal(.rowMax(m), sparseVector(c(5, 10, 15, 20), 1:4, 4))
    expect_equal(as.vector(.rowMax(m)), c(5, 10, 15, 20))
    expect_equal(as.vector(.rowMax(t(m))), 11:20)

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    expect_equal(as.vector(.rowMax(n, na.rm=FALSE)),
                           c(rep(NA, 3), 20))
    expect_equal(as.vector(.rowMax(n, na.rm=TRUE)), c(5, 10, 15, 20))
})

test_that(".rowMeansGroup", {
    r <- sparseMatrix(i=1:4,
                      j=c(1, 2, 1, 2),
                      x=c(3, 8, 13, 18))
    r2 <- sparseMatrix(i=1:10,
                       j=rep(1:2, each=5),
                       x=6:15)
    expect_error(.rowMeansGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(.rowMeansGroup(m, group=1:2))
    expect_equal(.rowMeansGroup(m, group=rep(1:2, each=5)), r)
    expect_equal(.rowMeansGroup(t(m), group=rep(1:2, 2)), r2)

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    r[cbind(1:3, c(1, 2, 1))] <- NA
    expect_equal(.rowMeansGroup(n, group=rep(1:2, each=5),
                                          na.rm=FALSE), r)
    r[cbind(1:3, c(1, 2, 1))] <- c(14/4, 25/3, 14)
    expect_equal(.rowMeansGroup(n, group=rep(1:2, each=5),
                                          na.rm=TRUE), r)
})

test_that(".rowSdsGroup", {
    .rowsds <- function(x, group, na.rm=TRUE) {
        l <- lapply(split(1L:ncol(x), group), function(i) {
            apply(x[, i, drop=FALSE], 1, sd, na.rm=na.rm)
        })
        do.call(cbind, l)
    }
    bm <- as.matrix(m)
    bm[bm == 0L] <- NA
    group <- rep(1:5, each=2)
    bmr <- .rowsds(bm, group)
    dimnames(bmr) <- NULL
    mr <- as.matrix(.rowSdsGroup(m, group))
    mr[mr == 0L] <- NA
    dimnames(mr) <- NULL
    expect_equal(bmr, mr)

    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    bn <- as.matrix(n)
    bn[bn == 0L] <- NA
    group <- rep(1:2, each=5)
    bnr <- .rowsds(bn, group, na.rm=FALSE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(.rowSdsGroup(n, group, na.rm=FALSE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)

    bnr <- .rowsds(bn, group, na.rm=TRUE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(.rowSdsGroup(n, group, na.rm=TRUE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)
})

test_that(".rowSumsGroup", {
    r <- sparseMatrix(i=rep(1:4, each=3),
                      j=rep(c(1:3, 3:5), 2),
                      x=c(3, 7, 5, 6, 15, 19, 23, 27, 15, 16, 35, 39))
    expect_error(.rowSumsGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(.rowSumsGroup(m, group=1:2))
    expect_equal(.rowSumsGroup(m, group=rep(1:5, each=2)), r)
    expect_equal(t(.rowSumsGroup(t(m), group=rep(1:2, 2))),
                 .colSumsGroup(m, group=rep(1:2, 2)))

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    r[cbind(c(1, 3, 2), c(1, 1, 4))] <- NA
    expect_equal(.rowSumsGroup(n, group=rep(1:5, each=2),
                                         na.rm=FALSE), r)
    r[cbind(c(1, 3, 2), c(1, 1, 4))] <- c(2, 0, 0)
    r <- drop0(r)
    expect_equal(.rowSumsGroup(n, group=rep(1:5, each=2),
                                         na.rm=TRUE), r)
})

test_that(".summary", {
    mr <- seq(1, 20, by=5)
    dr <- data.frame(Fragments=5, Total=c(15, 40, 65, 90),
                     Min=mr, Q1=mr + 1, Median=mr + 2, Mean=(mr + 2)/2,
                     Q3=mr + 3, Max=mr + 4)
    mc <- 1:10
    dc <- data.frame(Fragments=2, Total=mc * 2 + 10,
                     Min=mc, Q1=mc + 2.5, Median=mc + 5, Mean=(mc + 5)/2,
                     Q3=mc + 7.5, Max=mc + 10)
    expect_equal(.summary(m), dr)
    expect_equal(.summary(m, "row"), dr)
    expect_equal(.summary(m, "col"), dc)
})
