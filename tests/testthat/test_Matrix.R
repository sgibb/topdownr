context("Matrix")

test_that(".createMaskMatrix", {
  r1 <- sparseMatrix(i=1:10, j=rep(1:5, 2), x=1)
  r2 <- sparseMatrix(i=1:10, j=rep(1:2, each=5), x=1)
  expect_equal(topdown:::.createMaskMatrix(rep(1:5, 2)), r1)
  expect_equal(topdown:::.createMaskMatrix(rep(letters[1:5], 2)), r1)
  expect_equal(topdown:::.createMaskMatrix(rep(1:2, each=5)), r2)
})

test_that(".rowMeansGroup", {
  m <- sparseMatrix(i=rep(1:4, each=5),
                    j=rep(1:10, 2),
                    x=1:20)
  r <- sparseMatrix(i=c(1:4),
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

test_that(".colSumsGroup", {
  m <- sparseMatrix(i=rep(1:4, each=5),
                    j=rep(1:10, 2),
                    x=1:20)
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
