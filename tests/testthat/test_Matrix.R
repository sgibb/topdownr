context("Matrix")

test_that(".rowMeansGroup", {
  m <- sparseMatrix(i=rep(1:4, each=5),
                    j=rep(1:10, 2),
                    x=1:20)
  r <- sparseMatrix(i=c(1:4),
                    j=c(1, 2, 1, 2),
                    x=c(3, 8, 13, 18))
  r2 <- sparseMatrix(i=rep(1:10, 2),
                     j=rep(1:2, each=10),
                     x=c(seq(0.5, 10, by=0.5)))
  expect_error(topdown:::.rowMeansGroup(m, group=1:2))
  expect_equal(topdown:::.rowMeansGroup(m, group=rep(1:2, each=5)), r)
  expect_equal(topdown:::.rowMeansGroup(t(m), group=rep(1:2, each=2)), r2)
})
