context("Spectrum2")

test_that(".subsetSpectrum2", {
  s <- new("Spectrum2", mz=1:5, intensity=5:1,
           acquisitionNum=1L, fromFile=1L)
  r <- new("Spectrum2", mz=3:5, intensity=3:1,
           acquisitionNum=1L, fromFile=1L)
  e <- new("Spectrum2", fromFile=1L, acquisitionNum=1L)
  expect_error(topdown:::.subsetSpectrum2(s, 6))
  expect_equal(topdown:::.subsetSpectrum2(s, 3:5), r)
  expect_equal(topdown:::.subsetSpectrum2(s, integer()), e)
})
