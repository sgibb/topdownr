context("method-xml")

test_that("massLabel", {
  expect_equal(topdown:::massLabel(c(750, 1000.76), c(1, 245)),
               c(750.0001, 1000.8245))
  expect_equal(topdown:::massLabel(c(750, 1000.76), c(1, 245), divisor=1e5),
               c(750.00001, 1000.80245))
  expect_error(topdown:::massLabel(c(750, 1000.76), c(1, 245), divisor=1e3),
               "at least two digits more than")
})
