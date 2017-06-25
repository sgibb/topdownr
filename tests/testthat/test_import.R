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
                  TargetedMassList=paste0("(mz=933.100", 1:3, " z=2 name=)"),
                  stringsAsFactors=FALSE)
  write.csv(d, file=fn, row.names=FALSE)
  expect_message(e <- topdown:::.readExperimentCsv(fn, verbose=TRUE),
                 "Reading 3 experiment conditions from file")
  expect_equal(colnames(e),
               c("MSLevel", "TargetedMassList", "ConditionId", "Mz", "File"))
  expect_equal(e$MSLevel, rep(2, 2))
  expect_equal(e$ConditionId, 1:2)
  expect_equal(e$Mz, rep(933.1, 2))
  expect_equal(e$File, rep(basename(fn), 2))
  unlink(fn)
})

test_that(".readScanHeadsTable", {
  fn <- paste0(tempfile(), ".txt")
  d <- data.frame(MSOrder=c(1, 2, 2),
                  stringsAsFactors=FALSE)
  write.csv(d, file=fn, row.names=FALSE)
  expect_message(e <- topdown:::.readScanHeadsTable(fn, verbose=TRUE),
                 "Reading 3 header information from file")
  expect_equal(colnames(e),
               c("MSOrder", "File"))
  expect_equal(e$MSOrder, rep(2, 2))
  expect_equal(e$File, rep(basename(fn), 2))
  unlink(fn)
})
