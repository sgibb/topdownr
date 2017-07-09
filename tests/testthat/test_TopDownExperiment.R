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

## filtered TDE
f <- new.env()
f$F1.S1 <- e$F1.S1
f$F1.S2 <- new("Spectrum2", mz=5, intensity=9,
               acquisitionNum=2L, fromFile=1L)
f$F1.S3 <- new("Spectrum2", mz=c(3, 5), intensity=c(9, 5),
               acquisitionNum=3L, fromFile=1L)
atabf <- data.table(SpectrumId=paste0("F1.S", rep(1:3, c(3, 1, 2))),
                    FragmentId=c(1:3, 3, 2:3), MzId=c(1:3, 1, 1:2),
                    key=c("SpectrumId", "FragmentId", "MzId"))
tdf <- td
tdf@assayData <- f
tdf@assignmentTable <- atabf

f2 <- new.env()
f2$F1.S1 <- new("Spectrum2", mz=c(1, 3), intensity=c(5, 3),
                acquisitionNum=1L, fromFile=1L)
f2$F1.S2 <- new("Spectrum2", mz=7, intensity=7,
                acquisitionNum=2L, fromFile=1L)
f2$F1.S3 <- new("Spectrum2", mz=3, intensity=9,
                acquisitionNum=3L, fromFile=1L)
atabf2 <- data.table(SpectrumId=paste0("F1.S", rep(1:3, c(2, 1, 1))),
                     FragmentId=c(1:2, 4, 2), MzId=c(1:2, 1, 1),
                     key=c("SpectrumId", "FragmentId", "MzId"))
tdf2 <- td
tdf2@assayData <- f2
tdf2@assignmentTable <- atabf2

test_that("[", {
  e1 <- new.env()
  e1$F1.S1 <- e$F1.S1
  e12 <- new.env()
  e12$F1.S1 <- e$F1.S1
  e12$F1.S2 <- e$F1.S2
  td1 <- new("TopDownExperiment",
             assayData=e1,
             sequence="ACE",
             featureData=new("AnnotatedDataFrame", data=fd[1,,drop=FALSE]),
             phenoData=new("NAnnotatedDataFrame", data=pd),
             experimentData=expdata,
             processingData=new("MSnProcess", files="foo.mzML"),
             fragmentTable=ftab, assignmentTable=atab[1:3,])
  td12 <- new("TopDownExperiment",
              assayData=e12,
              sequence="ACE",
              featureData=new("AnnotatedDataFrame", data=fd[1:2,,drop=FALSE]),
              phenoData=new("NAnnotatedDataFrame", data=pd),
              experimentData=expdata,
              processingData=new("MSnProcess", files="foo.mzML"),
              fragmentTable=ftab, assignmentTable=atab[1:6,])

  expect_equal_TDE(td[1], td1)
  expect_equal_TDE(td["F1.S1"], td1)
  expect_equal_TDE(td[c(TRUE, FALSE, FALSE)], td1)
  expect_equal_TDE(td[1:2], td12)
  expect_equal_TDE(td[c("F1.S1", "F1.S2")], td12)
  expect_equal_TDE(td[c(TRUE, TRUE, FALSE)], td12)
  expect_equal_TDE(td[1, c("b", "c")], td1)
  expect_equal_TDE(td[1:2, c("b", "c")], td12)
  expect_equal_TDE(td[1:2, 1:2], td12)
  expect_equal_TDE(td[, "b"], tdf)
  expect_equal_TDE(td[, c("b1", "b2")], tdf)
  expect_equal_TDE(td[, c("b", "b2")], tdf)
  expect_equal_TDE(td[, 1], tdf2)
  expect_true(grepl("Subset [3;9] to [3;4]",
                    processingData(td[, 1])@processing, fixed=TRUE))
})

test_that(".matchFragments", {
  ## do nothing
  expect_equal(topdown:::.matchFragments(td, ftab, verbose=FALSE),
               list(assayData=assayData(td), assignmentTable=td@assignmentTable))

  eb <- new.env()
  eb$F1.S1 <- e$F1.S1
  eb$F1.S2 <- new("Spectrum2", mz=5, intensity=9,
                  acquisitionNum=2L, fromFile=1L)
  eb$F1.S3 <- new("Spectrum2", mz=c(3, 5), intensity=c(9, 5),
                  acquisitionNum=3L, fromFile=1L)
  ftabb <- data.table(mz=c(1, 3, 5),
                      ion=c("b1", "b1", "b2"),
                      type=c("b", "b", "b"),
                      pos=c(1, 1, 2),
                      z=1, FragmentId=1:3, key="FragmentId")
  atabb <- data.table(SpectrumId=paste0("F1.S", rep(1:3, c(3, 1, 2))),
                      FragmentId=c(1:3,
                                   3,
                                   2:3),
                      MzId=c(1:3,
                             1,
                             1:2), key=c("SpectrumId", "FragmentId", "MzId"))

  expect_equal(topdown:::.matchFragments(td, ftabb, verbose=FALSE),
               list(assayData=eb, assignmentTable=atabb))
})

test_that("assignmentTable", {
  expect_equal(topdown:::assignmentTable(td), atab)
})

test_that("dim", {
  expect_equal(dim(td), as.integer(c(length(td), nrow(atab))))
})

test_that(".filterFragmentId", {
  expect_error(topdown:::.filterFragmentId(td, 25), "empty object")
  expect_equal_TDE(topdown:::.filterFragmentId(td, 1:5), td)
  expect_equal_TDE(topdown:::.filterFragmentId(td, 1:3), tdf)
})

test_that(".filterFragmentIonOrType", {
  expect_error(topdown:::.filterFragmentIonOrType(td, c("b", "J", "D")),
               "Ion\\(s\\)/Type\\(s\\) .*J.*, .*D.* not found")
  expect_equal_TDE(topdown:::.filterFragmentIonOrType(td, paste0(c("b", "c"),
                                                                 c(1:2, 2:1))),
                   td)
  expect_equal_TDE(topdown:::.filterFragmentIonOrType(td, paste0("b", 1:2)), tdf)
  expect_equal_TDE(topdown:::.filterFragmentIonOrType(td, "b"), tdf)
})

test_that(".filterFragmentPos", {
  expect_error(topdown:::.filterFragmentPos(td, c(100, 1000)),
               "Position\\(s\\) .*100.*, .*1000.* not found")
  expect_equal_TDE(topdown:::.filterFragmentPos(td, 1:2), td)
  expect_equal_TDE(topdown:::.filterFragmentPos(td, 1), tdf2)
})

test_that("fragmentTable", {
  expect_equal(topdown:::fragmentTable(td), ftab)
})

test_that(".fragmentId", {
  expect_equal(topdown:::.fragmentId(td), ftab$FragmentId)
})

test_that(".fragmentPos", {
  expect_equal(topdown:::.fragmentPos(td), c(1, 1, 2, 1, 2))
})

test_that(".fragmentTypes", {
  expect_equal(topdown:::.fragmentTypes(td), c("b", "b", "b", "c", "c"))
})

test_that(".logmsg", {
  expect_equal(topdown:::.logmsg(td, "foo",
                                 date=FALSE)@processingData@processing[1], "foo")
  expect_equal(topdown:::.logmsg(td, "foo")@processingData@processing[1],
               paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] foo"))
})

test_that(".validateTopDownExperiment", {
  expect_true(topdown:::.validateTopDownExperiment(td))
  expect_true(validObject(td))
  td@assignmentTable$SpectrumId <- paste0("F2.S", 1:9)
  expect_error(validObject(td),
               "IDs in assignment table don't match feature names.")
  td@assignmentTable <- atab
  td@assignmentTable$MzId <- 1:9
  expect_error(validObject(td),
               "Mismatch in spectra and assignment table's peak indices.")
  td@assignmentTable <- atab
  td@assignmentTable$FragmentId <- 1:9
  expect_error(validObject(td),
               "Mismatch in fragment and assignment table's fragment indices.")
})
