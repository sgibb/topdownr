context("utils")

test_that("cat0", {
  expect_output(topdown:::cat0("foo", "bar"), "foobar")
})

test_that(".massLabel", {
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245)),
               c(750.0001, 1000.8245))
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e5),
               c(750.00001, 1000.80245))
  expect_error(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e3),
               "at least two digits more than")
})

test_that(".massLabelToId", {
  expect_equal(topdown:::.massLabelToId(c("750.0001", "1000.8245")), c(1, 245))
  expect_equal(topdown:::.massLabelToId(c(750.0001, 1000.8245)), c(1, 245))
  expect_equal(topdown:::.massLabelToId(topdown:::.massLabel(c(750, 1000.76),
                                                             c(1, 245))),
                                        c(1, 245))
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

test_that(".targetMassListToMz", {
  expect_error(topdown:::.targetMassListToMz(1:3))
  expect_equal(topdown:::.targetMassListToMz(c("(mz=1000.12 z=2 name=foo)",
                                               "(mz=933.99 z=3 name=)")),
               c(1000.1, 933.9))
})
