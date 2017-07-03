context("TopDownExperiment")

expect_equal_TDE <- function(object, expected, ..., proc=FALSE, index=FALSE,
                             info=NULL, label=NULL) {
  if (!proc) {
    object@processingData <-
      expected@processingData <-
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
atab <- data.table(SpectrumId=paste0("F1.S", rep(1:3, each=3)),
                   FragmentId=c(1:3,
                                3:5,
                                2:3, 5),
                   MzId=c(1:3,
                          1:3,
                          1:3), key=c("SpectrumId", "FragmentId", "MzId"))
td <- new("TopDownExperiment",
          assayData=e,
          sequence="ACE",
          featureData=new("AnnotatedDataFrame", data=fd),
          phenoData=new("NAnnotatedDataFrame", data=pd),
          experimentData=expdata,
          processingData=new("MSnProcess", files="foo.mzML"),
          fragmentTable=ftab, assignmentTable=atab)

test_that("[", {
  e1 <- new.env()
  e1$F1.S1 <- e$F1.S1
  td1 <- new("TopDownExperiment",
             assayData=e1,
             sequence="ACE",
             featureData=new("AnnotatedDataFrame", data=fd[1,,drop=FALSE]),
             phenoData=new("NAnnotatedDataFrame", data=pd[1,,drop=FALSE]),
             experimentData=expdata,
             processingData=new("MSnProcess", files="foo.mzML"),
             fragmentTable=ftab, assignmentTable=atab[1:3,])
  expect_equal_TDE(td[1], td1)
})

test_that("assignmentTable", {
  expect_equal(topdown:::assignmentTable(td), atab)
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
