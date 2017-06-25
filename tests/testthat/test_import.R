context("import")

test_that(".fileExists", {
  expect_error(topdown:::.fileExists("foo.bar"))
  fn <- tempfile()
  file.create(fn)
  expect_true(topdown:::.fileExists(fn))
  unlink(fn)
})

test_that(".listTopDownFiles", {
  fns <- tempfile(pattern=c("fileA_", "fileB_"))
  fns <- paste(rep(fns, each=3), c("experiments.csv", "mzML", "txt"), sep=".")
  r <- split(fns, c("csv", "mzML", "txt"))
  file.create(fns)
  expect_equal(topdown:::.listTopDownFiles(tempdir()), r)
  expect_equal(topdown:::.listTopDownFiles(tempdir(), pattern="^fileA_.*"),
               lapply(r, "[", 1L))
  unlink(fns)
})

test_that(".readExperimentCsv", {
  fn <- paste0(tempfile(), ".experiments.csv")
  d <- data.frame(MSLevel=c(1, 2, 2),
                  NaColumn=NA,
                  TargetedMassList=paste0("(mz=933.100", 1:3, " z=2 name=)"),
                  stringsAsFactors=FALSE)
  write.csv(d, file=fn, row.names=FALSE)
  expect_message(e <- topdown:::.readExperimentCsv(fn, verbose=TRUE),
                 "Reading 3 experiment conditions from file")
  expect_equal(colnames(e),
               c("MSLevel", "NaColumn", "TargetedMassList", "ConditionId",
                 "Mz", "File"))
  expect_equal(e$MSLevel, rep(2, 2))
  expect_equal(e$NaColumn, rep(0, 2))
  expect_equal(e$ConditionId, 1:2)
  expect_equal(e$Mz, rep(933.1, 2))
  expect_equal(e$File, rep(basename(fn), 2))
  unlink(fn)
})

test_that(".readScanHeadsTable", {
  fn <- paste0(tempfile(), ".txt")
  d <- data.frame(MSOrder=c(1, 2, 2, 2),
                  FilterString=c("ms2 100.0001@etd", "ms2 100.0001@hcd",
                                 "ms2 100.0007@cid", "ms2 100.0009@hcd"),
                  Activation1=c("ETD", "ETD", "HCD", "CID"),
                  Activation2=c(NA, "HCD", "CID", "HCD"),
                  Energy1=c(10, 50, NA, 20),
                  Energy2=c(NA, 30, 20, 10),
                  stringsAsFactors=FALSE)
  write.csv(d, file=fn, row.names=FALSE)
  expect_message(h <- topdown:::.readScanHeadsTable(fn, verbose=TRUE),
                 "Reading 4 header information from file")
  expect_equal(colnames(h),
               c("MSOrder", "FilterString", "Activation1", "Activation2",
                 "Energy1", "Energy2", "ETDActivation", "CIDActivation",
                 "HCDActivation", "ConditionId", "File"))
  expect_equal(h$MSOrder, rep(2, 3))
  expect_equal(h$ETDActivation, c(50, 0, 0))
  expect_equal(h$CIDActivation, c(0, 20, 20))
  expect_equal(h$HCDActivation, c(30, 0, 10))
  expect_equal(h$ConditionId, c(1, 7, 9))
  expect_equal(h$File, rep(basename(fn), 3))
  unlink(fn)
})
