context("FragmentViews")

fv <- FragmentViews("ACE", start=3:1, end=3, names=paste0("x", 1:3),
                    mass=c(174.039696, 334.070346, 405.107456),
                    type=rep("x", 3), z=rep(1, 3))

test_that("constructor", {
  expect_equal(topdown:::.calculateFragments("ACE", type="x",
                                             neutralLoss=list(water=NULL)), fv)
})

test_that("show", {
  expect_output(show(fv),
                paste(c("FragmentViews on a 3-letter sequence:",
                        "  ACE",
                        "Views:",
                        "    start end width   mass type z *",
                        "\\[1\\]     3   3     1 174\\.04 x    1 \\[E\\] *",
                        "\\[2\\]     2   3     2 334\\.07 x    1 \\[CE\\] *",
                        "\\[3\\]     1   3     3 405\\.11 x    1 \\[ACE\\]"),
                      collapse="\n"))
})

test_that("validity", {
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
