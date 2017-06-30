context("MSnExp")

expect_equal_MSnExp <- function(object, expected, ..., proc=FALSE,
                             info=NULL, label=NULL) {
  if (!proc) {
    object@processingData <- expected@processingData <- new("MSnProcess")
  }
  expect_equal(object, expected, ..., info=info, label=label)
}

test_that(".subsetMSnExpSpectra", {
  e <- new.env()
  e$F1.S1 <- new("Spectrum2", mz=1:5, intensity=5:1,
                 acquisitionNum=1L, fromFile=1L)
  e$F1.S2 <- new("Spectrum2", mz=4:9, intensity=9:4,
                 acquisitionNum=2L, fromFile=1L)
  e$F1.S3 <- new("Spectrum2", mz=2:6, intensity=6:2,
                 acquisitionNum=3L, fromFile=1L)
  fd <- data.frame(fileIdx=rep(1, 3),
                   Scan=1:3,
                   File="foo",
                   spectrum=1:3,
                   row.names=paste0("F1.S", 1:3))
  pd <- data.frame(sampleNames = "foo.mzML",
                   row.names="foo.mzML")
  expdata <- new("MIAPE",
                 instrumentManufacturer = "Thermo",
                 instrumentModel = "Orbitrap Fusion Luminos",
                 ionSource = "nanoelectrospray",
                 analyser = "orbitrap",
                 detectorType = "Unknown")
  mse <- new("MSnExp", assayData=e,
             featureData=new("AnnotatedDataFrame", data=fd),
             phenoData=new("NAnnotatedDataFrame", data=pd),
             experimentData=expdata,
             processingData=new("MSnProcess", files="foo.mzML"))

  expect_error(topdown:::.subsetMSnExpSpectra(mse, 6, list(1, 1)))
  expect_equal_MSnExp(topdown:::.subsetMSnExpSpectra(mse, rep(1:3, c(5, 6, 5)),
                                                     c(1:5, 1:6, 1:5)), mse)
})
