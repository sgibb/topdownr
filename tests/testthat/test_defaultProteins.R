context("defaultProteins")

test_that("defaultProteins", {
  nm <- c("GST", "myoglobin", "h2a", "h2b", "h3_1", "h3_3", "h4",
          "c345c_c3", "h33tail")
  expect_equal(names(defaultProteins), nm)
  expect_equal(defaultProteins$h2a,
               matrix(c(609.21, 700.45, 823.95, 10, 10, 10),
                      ncol=2, dimnames=list(c(), c("mass", "z"))))
})

test_that(".mzMatrix", {
  m <- matrix(c(100, 200, 300, 10, 10, 10),
              ncol=2, dimnames=list(c(), c("mass", "z")))
  expect_equal(topdown:::.mzMatrix(c(100, 200, 300)), m)
  expect_equal(topdown:::.mzMatrix(100), m[1,, drop=FALSE])
})

