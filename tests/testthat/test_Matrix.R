context("Matrix")

m <- sparseMatrix(i=rep(1:4, each=5), j=rep(1:10, 2), x=1:20)

test_that(".bestNcbCoverageCombination", {
    m1 <- sparseMatrix(i=c(1, 2, 2, 3, 5, 5, 6, 9, 8, 9, 8),
                       j=c(1, 5, 3, 4, 4, 5, 4, 4, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 2))
    expect_equal(topdown:::.bestNcbCoverageCombination(m1),
                 cbind(index=c(5:4, 1, 3), fragments=c(5:4, 1, 1)))
    expect_equal(topdown:::.bestNcbCoverageCombination(m1, intensity=1:5),
                 cbind(index=c(5:3, 1), fragments=c(5:4, 1, 1)))
    expect_equal(topdown:::.bestNcbCoverageCombination(m1, minN=3),
                 cbind(index=c(5:4), fragments=c(5, 4)))
    expect_equal(topdown:::.bestNcbCoverageCombination(m1, n=3),
                 cbind(index=c(5:4, 1), fragments=c(5:4, 1)))
    expect_equal(topdown:::.bestNcbCoverageCombination(m1, minN=5),
                 cbind(index=5, fragments=5))
})

test_that(".col", {
    expect_error(topdown:::.col(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.col(m), rep(1:10, each=2))
})

test_that(".colCounts", {
    expect_error(topdown:::.colCounts(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.colCounts(m), rep(2, 10))
})

test_that(".colSumsGroup", {
    r <- sparseMatrix(i=rep(1:2, each=5),
                      j=1:10,
                      x=seq(12, 30, by=2))
    r2 <- sparseMatrix(i=rep(1:2, 2),
                       j=1:4,
                       x=c(15, 40, 65, 90))
    expect_error(topdown:::.colSumsGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(topdown:::.colSumsGroup(m, group=1:2))
    expect_equal(topdown:::.colSumsGroup(m, group=rep(1:2, 2)), r)
    expect_equal(topdown:::.colSumsGroup(t(m), group=rep(1:2, each=5)), r2)
})

test_that(".countFragments", {
    expect_error(topdown:::.countFragments(matrix(1:10, ncol=2)))
    expect_error(topdown:::.countFragments(m))
    expect_equal(topdown:::.countFragments(drop0(m %% 4)), rep(c(3, 1), 5))
})

test_that(".createMaskMatrix", {
    r1 <- sparseMatrix(i=1:10, j=rep(1:5, 2), x=1)
    r2 <- sparseMatrix(i=1:10, j=rep(1:2, each=5), x=1)
    expect_equal(topdown:::.createMaskMatrix(rep(1:5, 2)), r1)
    expect_equal(topdown:::.createMaskMatrix(rep(letters[1:5], 2)), r1)
    expect_equal(topdown:::.createMaskMatrix(rep(1:2, each=5)), r2)
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
    expect_error(topdown:::.cumComb(1:10))
    expect_equal(topdown:::.cumComb(a), r)
    expect_equal(topdown:::.cumComb(a2), a2)
})

test_that(".dgcMatrix2data.frame", {
    d <- data.frame(row=c(rep(c(1, 3), 5), rep(c(2, 4), 5)),
                    col=rep(1:10, each=2), x=rep(1:10, each=2) + c(0, 10),
                    stringsAsFactors=FALSE)
    expect_error(topdown:::.dgcMatrix2data.frame(matrix(1:10, nrow=2)))
    expect_equal(topdown:::.dgcMatrix2data.frame(m), d)
})

test_that(".drop0rowLe/Lt", {
    r1 <- sparseMatrix(i=rep(1:4, each=3), j=rep(c(3:5, 8:10), 2),
                       x=c(3:5, 8:10, 13:15, 18:20))
    r2 <- sparseMatrix(i=1:10, j=rep(3:4, each=5), x=11:20)
    expect_error(topdown:::.drop0rowLe(matrix(1:10, ncol=2), 1))
    expect_error(topdown:::.drop0rowLt(matrix(1:10, ncol=2), 1))
    expect_error(topdown:::.drop0rowLe(m, 1))
    expect_error(topdown:::.drop0rowLt(m, 1))
    expect_equal(topdown:::.drop0rowLe(m, c(2, 7, 12, 17)), r1)
    expect_equal(topdown:::.drop0rowLt(m, c(3, 8, 13, 18)), r1)
    expect_equal(topdown:::.drop0rowLe(t(m), 1:10), r2)
    expect_equal(topdown:::.drop0rowLt(t(m), 2:11), r2)
})

test_that(".drop0rowReplicates", {
    group <- rep(1:5, each=2)
    expect_error(topdown:::.drop0rowReplicates(matrix(1:10, ncol=2), 1))
    expect_error(topdown:::.drop0rowReplicates(m, 1))
    expect_error(topdown:::.drop0rowReplicates(m, group, 1.5))
    expect_error(topdown:::.drop0rowReplicates(m, group, 1L:2L))
    expect_equal(topdown:::.drop0rowReplicates(m, group, minN=1L), m)
    expect_equal(topdown:::.drop0rowReplicates(m, group, minN=2L),
                 sparseMatrix(i=rep(1:4, each=4),
                              j=rep(c(1:4, 7:10), 2),
                              x=c(1:4, 7:10, 11:14, 17:20)))
    expect_equal(topdown:::.drop0rowReplicates(m, group, minN=3L),
                 Matrix(0L, nrow=4, ncol=10, sparse=TRUE))

    expect_equal(topdown:::.drop0rowReplicates(
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
    expect_error(topdown:::.dropNA(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.dropNA(m), m)
    expect_equal(topdown:::.dropNA(n), r)
})

test_that(".highestNcbCoverage", {
    m1 <- sparseMatrix(i=c(1:10, 2), j=c(1, 5, 3, 4, 4, 5, 4, 5, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 1))
    expect_error(topdown:::.highestNcbCoverage(as.matrix(1:10)))
    expect_error(topdown:::.highestNcbCoverage(m))
    expect_equal(topdown:::.highestNcbCoverage(m1, intensity=1:5),
                 c(index=5, fragments=6))
    expect_equal(topdown:::.highestNcbCoverage(t(m1)),
                 c(index=2, fragments=2))
})

test_that(".normaliseCols", {
    expect_error(topdown:::.normaliseCols(matrix(1:10, nrow=2)))
    expect_error(topdown:::.normaliseCols(m, "A"))
    expect_error(topdown:::.normaliseCols(m, 1:2))
    expect_equal(topdown:::.normaliseCols(m),
                 as((scale(m, center=FALSE, scale=11:20)),
                    "dgCMatrix"))
    expect_equal(topdown:::.normaliseCols(m, 11:20),
                 as((scale(m, center=FALSE, scale=11:20)),
                    "dgCMatrix"))
})

test_that(".normaliseRows", {
    expect_error(topdown:::.normaliseRows(matrix(1:10, nrow=2)))
    expect_error(topdown:::.normaliseRows(m, "A"))
    expect_error(topdown:::.normaliseRows(m, 1:2))
    expect_equal(topdown:::.normaliseRows(m),
                 as((t(scale(t(m), center=FALSE, scale=c(5, 10, 15, 20)))),
                    "dgCMatrix"))
    expect_equal(topdown:::.normaliseRows(m, 1:4),
                 as((t(scale(t(m), center=FALSE, scale=1:4))),
                    "dgCMatrix"))
})

test_that(".removeNcbCombinations", {
    m1 <- sparseMatrix(i=c(1, 2, 2, 3, 5, 5, 6, 9, 8, 9, 8),
                       j=c(1, 5, 3, 4, 4, 5, 4, 4, 5, 5, 1),
                       x=c(1, 1, 3, 1, 1, 2, 2, 2, 3, 1, 2))
    r1 <- sparseMatrix(i=c(1, 2, 3, 5, 6, 9),
                       j=c(1, 3, 4, 4, 4, 4),
                       x=c(1, 2, 1, 1, 2, 2), dims=c(9, 5))
    r2 <- sparseMatrix(i=c(1, 2, 2, 5, 8, 9, 8),
                       j=c(1, 5, 3, 5, 5, 5, 1),
                       x=c(1, 1, 3, 2, 3, 1, 2))
    expect_error(topdown:::.removeNcbCombinations(matrix(1:10, ncol=2), 2))
    expect_error(topdown:::.removeNcbCombinations(m, 2))
    expect_equal(topdown:::.removeNcbCombinations(m1, 5), r1)
    expect_equal(topdown:::.removeNcbCombinations(m1, 4), r2)
})

test_that(".row", {
    expect_error(topdown:::.row(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.row(m), c(rep(c(1, 3), 5), rep(c(2, 4), 5)))
})

test_that(".rowCounts", {
    expect_error(topdown:::.rowCounts(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.rowCounts(m), rep(5, 4))
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
    mr <- as.matrix(topdown:::.rowCvsGroup(m, group))
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
    nr <- as.matrix(topdown:::.rowCvsGroup(n, group, na.rm=FALSE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)

    bnr <- .rowcvs(bn, group, na.rm=TRUE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(topdown:::.rowCvsGroup(n, group, na.rm=TRUE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)
})

test_that(".rowMax", {
    expect_error(topdown:::.rowMax(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.rowMax(m), sparseVector(c(5, 10, 15, 20), 1:4, 4))
    expect_equal(as.vector(topdown:::.rowMax(m)), c(5, 10, 15, 20))
    expect_equal(as.vector(topdown:::.rowMax(t(m))), 11:20)

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    expect_equal(as.vector(topdown:::.rowMax(n, na.rm=FALSE)),
                           c(rep(NA, 3), 20))
    expect_equal(as.vector(topdown:::.rowMax(n, na.rm=TRUE)), c(5, 10, 15, 20))
})

test_that(".rowMeansGroup", {
    r <- sparseMatrix(i=1:4,
                      j=c(1, 2, 1, 2),
                      x=c(3, 8, 13, 18))
    r2 <- sparseMatrix(i=1:10,
                       j=rep(1:2, each=5),
                       x=6:15)
    expect_error(topdown:::.rowMeansGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(topdown:::.rowMeansGroup(m, group=1:2))
    expect_equal(topdown:::.rowMeansGroup(m, group=rep(1:2, each=5)), r)
    expect_equal(topdown:::.rowMeansGroup(t(m), group=rep(1:2, 2)), r2)

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    r[cbind(1:3, c(1, 2, 1))] <- NA
    expect_equal(topdown:::.rowMeansGroup(n, group=rep(1:2, each=5),
                                          na.rm=FALSE), r)
    r[cbind(1:3, c(1, 2, 1))] <- c(14/4, 25/3, 14)
    expect_equal(topdown:::.rowMeansGroup(n, group=rep(1:2, each=5),
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
    mr <- as.matrix(topdown:::.rowSdsGroup(m, group))
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
    nr <- as.matrix(topdown:::.rowSdsGroup(n, group, na.rm=FALSE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)

    bnr <- .rowsds(bn, group, na.rm=TRUE)
    dimnames(bnr) <- NULL
    nr <- as.matrix(topdown:::.rowSdsGroup(n, group, na.rm=TRUE))
    nr[nr == 0L] <- NA
    dimnames(nr) <- NULL
    expect_equal(bnr, nr)
})

test_that(".rowSumsGroup", {
    r <- sparseMatrix(i=rep(1:4, each=3),
                      j=rep(c(1:3, 3:5), 2),
                      x=c(3, 7, 5, 6, 15, 19, 23, 27, 15, 16, 35, 39))
    expect_error(topdown:::.rowSumsGroup(matrix(1:10, ncol=2), group=1:2))
    expect_error(topdown:::.rowSumsGroup(m, group=1:2))
    expect_equal(topdown:::.rowSumsGroup(m, group=rep(1:5, each=2)), r)
    expect_equal(t(topdown:::.rowSumsGroup(t(m), group=rep(1:2, 2))),
                 topdown:::.colSumsGroup(m, group=rep(1:2, 2)))

    ## na.rm
    n <- m
    n[cbind(c(1, 3, 3, 2, 2), c(1, 1, 2, 7, 8))] <- NA
    r[cbind(c(1, 3, 2), c(1, 1, 4))] <- NA
    expect_equal(topdown:::.rowSumsGroup(n, group=rep(1:5, each=2),
                                         na.rm=FALSE), r)
    r[cbind(c(1, 3, 2), c(1, 1, 4))] <- c(2, 0, 0)
    r <- drop0(r)
    expect_equal(topdown:::.rowSumsGroup(n, group=rep(1:5, each=2),
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
    expect_equal(topdown:::.summary(m), dr)
    expect_equal(topdown:::.summary(m, "row"), dr)
    expect_equal(topdown:::.summary(m, "col"), dc)
})
