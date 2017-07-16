context("fragments")

test_that(".matchFragments", {
  expect_equal(topdown:::.matchFragments(mz=integer(), fmass=1:3), integer())
  expect_equal(topdown:::.matchFragments(c(1, 99, 101), fmass=c(1.1, 100),
                                         tolerance=0.2),
               as.integer(c(1, 2, NA)))
  expect_equal(topdown:::.matchFragments(c(1, 98, 101), fmass=c(1.1, 100),
                                         tolerance=0.2),
               as.integer(c(1, NA, 2)))
})
