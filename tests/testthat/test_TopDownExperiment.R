context("TopDownExperiment")

expect_equal_TDE <- function(object, expected, ..., proc=FALSE, index=FALSE,
                             info=NULL, label=NULL) {
  if (!proc) {
    object@msnExp@processingData <-
      expected@msnExp@processingData <-
        new("MSnProcess")
  }
  if (!index) {
    attr(object@assignmentTable, "index") <- NULL
    attr(expected@assignmentTable, "index") <- NULL
  }
  expect_equal(object, expected, ..., info=info, label=label)
}

## basic TDE
e <- new.env()
e$F1.S1 <- new("Spectrum2", mz=c(1, 3, 5), intensity=c(5, 3, 1),
               acquisitionNum=1L, fromFile=1L)
e$F1.S2 <- new("Spectrum2", mz=c(5, 7, 9), intensity=c(9, 7, 5),
               acquisitionNum=2L, fromFile=1L)
e$F1.S3 <- new("Spectrum2", mz=c(3, 5, 9), intensity=c(9, 5, 3),
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
ftab <- data.table(mz=c(1, 3, 5, 7, 9),
                   ion=c("b1", "b1", "b2", "c1", "c2"),
                   type=c("b", "b", "b", "c", "c"),
                   pos=c(1, 1, 2, 1, 2),
                   z=1, FragmentId=1:5, key="FragmentId")
atab <- data.table(SpectrumId=rep(1:3, each=3),
                   FragmentId=c(1:3,
                                3:5,
                                2:3, 5),
                   MzId=c(1:3,
                          1:3,
                          1:3), key=c("SpectrumId", "FragmentId", "MzId"))
mse <- new("MSnExp", assayData=e,
           featureData=new("AnnotatedDataFrame", data=fd),
           phenoData=new("NAnnotatedDataFrame", data=pd),
           experimentData=expdata,
           processingData=new("MSnProcess", files="foo.mzML"))
td <- new("TopDownExperiment",
          sequence="ACE", msnExp=mse,
          fragmentTable=ftab, assignmentTable=atab)

## filtered TDE
f <- new.env()
f$F1.S1 <- e$F1.S1
f$F1.S2 <- new("Spectrum2", mz=5, intensity=9,
               acquisitionNum=2L, fromFile=1L)
f$F1.S3 <- new("Spectrum2", mz=c(3, 5), intensity=c(9, 5),
               acquisitionNum=3L, fromFile=1L)
atabf <- data.table(SpectrumId=c(1, 1, 1, 2, 3, 3),
                    FragmentId=c(1:3,
                                 3,
                                 2:3),
                    MzId=c(1:3,
                           1,
                           1:2), key=c("SpectrumId", "FragmentId", "MzId"))
msf <- mse
msf@assayData <- f
tdf <- new("TopDownExperiment",
           sequence="ACE", msnExp=msf,
           fragmentTable=ftab, assignmentTable=atabf)

test_that(".filterFragmentType", {
  expect_error(topdown:::.filterFragmentType(td, c("b", "J", "D")),
               "Type .*J.*, .*D.* is not valid")
})

test_that("assignmentTable", {
  expect_equal(topdown:::assignmentTable(td), atab)
})

test_that(".filterFragmentId", {
  expect_equal_TDE(topdown:::.filterFragmentId(td, 1:5), td)
  expect_equal_TDE(topdown:::.filterFragmentId(td, 1:3), tdf)
})

test_that(".filterFragmentType", {
  expect_equal_TDE(topdown:::.filterFragmentType(td, c("b", "c")), td)
  expect_equal_TDE(topdown:::.filterFragmentType(td, "b"), tdf)
})

test_that("fragmentTable", {
  expect_equal(topdown:::fragmentTable(td), ftab)
})

test_that(".fragmentId", {
  expect_equal(topdown:::.fragmentId(td), ftab$FragmentId)
})

test_that(".fragmentTypes", {
  expect_equal(topdown:::.fragmentTypes(td), c("b", "b", "b", "c", "c"))
})

test_that("msnExp", {
  expect_equal(msnExp(td), mse)
})
