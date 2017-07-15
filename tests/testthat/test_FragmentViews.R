context("FragmentViews")

test_that("constructor", {
  fv <- FragmentViews("ACE", start=3:1, end=3, names=paste0("x", 1:3),
                      mass=c(174.039696, 334.070346, 405.107456),
                      type=rep("x", 3), z=rep(1, 3))
  rv <- topdown:::.calculateFragments("ACE", type="x",
                                      neutralLoss=list(water=NULL))

})

test_that("validity", {
  fv <- FragmentViews("ACE", start=3:1, end=3, names=paste0("x", 1:3),
                      mass=as.double(1:3), type=rep("x", 3), z=rep(1, 3))
  expect_true(validObject(fv))
  d <- elementMetadata(fv)
  elementMetadata(fv) <- NULL
  expect_error(validObject(fv), "'mass', 'type', 'z' are missing")
  elementMetadata(fv) <- d
  elementMetadata(fv)$z <- NULL
  expect_error(validObject(fv), "'z' is missing")
  elementMetadata(fv) <- d
  elementMetadata(fv)$mass <- 1L:3L
  expect_error(validObject(fv), "'mass' has to be of type double")
  elementMetadata(fv) <- d
  elementMetadata(fv)$mass <- as.double(3:1)
  expect_error(validObject(fv), "'mass' has to be sorted")
})
