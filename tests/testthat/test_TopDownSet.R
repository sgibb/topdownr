context("TopDownSet")

test_that(".isTopDownSet", {
  expect_true(topdown:::.isTopDownSet(new("TopDownSet")))
  expect_error(topdown:::.isTopDownSet(1L), "has to be an 'TopDownSet' object")
})
