context("Spectrum2")

test_that(".aggregateSpectra", {
  mz <- list(A=1:5, B=1:3, C=3:4)
  int <- list(A=5:1, B=3:1, C=4:3)
  r <- new("Spectrum2", mz=1:5, intensity=c(4, 3, 8/3, 5/2, 1))

  ftab <- data.table(mz=1:7,
                     ion=c("a1", "a2", "b1", "b2", "c1", "c2", "c3"),
                     type=c("a", "a", "b", "b", "c", "c", "c"),
                     pos=c(1, 2, 1, 2, 1, 2, 3),
                     z=1, FragmentId=1:7, key="FragmentId")
  atab <- data.table(SpectrumId=rep(LETTERS[1:3], c(5, 3, 2)),
                     FragmentId=c(1:5,
                                  1:3,
                                  3:4),
                     MzId=c(1:5,
                            1:3,
                            1:2), key=c("SpectrumId", "FragmentId", "MzId"))

  expect_error(topdown:::.aggregateSpectra(unname(mz), int, atab, ftab))
  expect_error(topdown:::.aggregateSpectra(mz, int, atab[1:3,], ftab))
  expect_error(topdown:::.aggregateSpectra(mz[1], int, atab, ftab))
  expect_error(topdown:::.aggregateSpectra(mz, list(A=3:4, B=1:3, C=1:5), atab, ftab))
  expect_equal(topdown:::.aggregateSpectra(mz, int, atab, ftab), r)
})

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
