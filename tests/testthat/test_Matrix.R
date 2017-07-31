context("Matrix")

m <- sparseMatrix(i=rep(1:4, each=5), j=rep(1:10, 2), x=1:20)

test_that(".createMaskMatrix", {
    r1 <- sparseMatrix(i=1:10, j=rep(1:5, 2), x=1)
    r2 <- sparseMatrix(i=1:10, j=rep(1:2, each=5), x=1)
    expect_equal(topdown:::.createMaskMatrix(rep(1:5, 2)), r1)
    expect_equal(topdown:::.createMaskMatrix(rep(letters[1:5], 2)), r1)
    expect_equal(topdown:::.createMaskMatrix(rep(1:2, each=5)), r2)
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
                 sparseMatrix(i=rep(1:4, c(4, 2, 4, 2)),
                              j=rep(c(1:4, 9:10), 2),
                              x=c(1:4, 9:10, 11:14, 19:20)))
    expect_equal(topdown:::.drop0rowReplicates(m, group, minN=3L),
                 Matrix(0L, nrow=4, ncol=10, sparse=TRUE))
})

test_that(".m2rect", {
    m <- sparseMatrix(i=rep(1:2, each=5), j=1:10, x=1:10)
    r <- cbind(xleft=rep(0:1, each=5), ybottom=0:9,
               xright=rep(1:2, each=5), ytop=1:10,
               col=1:10)
    expect_error(topdown:::.m2rect(matrix(1:10, nrow=2)))
    expect_equal(topdown:::.m2rect(m), r)
})

test_that(".normaliseRows", {
    expect_error(topdown:::.normaliseRows(matrix(1:10, nrow=2)))
    expect_equal(topdown:::.normaliseRows(m),
                 as((t(scale(t(m), center=FALSE, scale=c(5, 10, 15, 20)))),
                    "dgCMatrix"))
})

test_that(".rowMax", {
    expect_error(topdown:::.rowMax(matrix(1:10, ncol=2)))
    expect_equal(topdown:::.rowMax(m), sparseVector(c(5, 10, 15, 20), 1:4, 4))
    expect_equal(as.vector(topdown:::.rowMax(m)), c(5, 10, 15, 20))
    expect_equal(as.vector(topdown:::.rowMax(t(m))), 11:20)
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

})
