context("utils")

test_that("cat0", {
  expect_output(topdown:::cat0("foo", "bar"), "foobar")
})

test_that(".msg", {
  expect_message(topdown:::.msg(TRUE, "foobar"), "foobar")
  expect_message(topdown:::.msg(TRUE, "foo", "bar"), "foobar")
  expect_silent(topdown:::.msg(FALSE, "foobar"))
})

test_that(".nrows", {
  expect_error(topdown:::.nrows(matrix(nrow=2, ncol=2)))
  expect_equal(topdown:::.nrows(list(matrix(nrow=2, ncol=2),
                                     matrix(nrow=3, ncol=2))), 2:3)
})

test_that(".targetMassList2Mz", {
  expect_error(topdown:::.targetMassList2Mz(1:3))
  expect_equal(topdown:::.targetMassList2Mz(c("(mz=1000.12 z=2 name=foo)",
                                              "(mz=933.99 z=3 name=)")),
               c(1000.1, 933.9))
})
