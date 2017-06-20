context("utils")

test_that("cat0", {
  expect_output(topdown:::cat0("foo", "bar"), "foobar")
})
