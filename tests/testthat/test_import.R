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
  fasta <- paste(tempfile(pattern=c("file")), "fasta", sep=".")
  fns <- paste(rep(fns, each=3), c("experiments.csv", "mzML", "txt"), sep=".")
  r <- split(fns, c("csv", "mzML", "txt"))
  r$fasta <- fasta
  r <- r[order(names(r))]
  file.create(c(fasta, fns))
  expect_equal(topdown:::.listTopDownFiles(tempdir()), r)
  expect_equal(topdown:::.listTopDownFiles(tempdir(), pattern="^fileA_.*"),
               lapply(r[names(r) != "fasta"], "[", 1L))
  unlink(c(fasta, fns))
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
               c("MSLevel", "NaColumn", "TargetedMassList", "Condition",
                 "Mz", "File"))
  expect_equal(e$MSLevel, rep(2, 2))
  expect_equal(e$NaColumn, rep(0, 2))
  expect_equal(e$Condition, 1:2)
  expect_equal(e$Mz, rep(933.1, 2))
  expect_equal(e$File, rep(gsub("\\.experiments.csv", "", basename(fn)), 2))
  unlink(fn)
})

test_that(".readFasta", {
  fn <- paste0(tempfile(), ".fasta")
  writeLines(c("> FOOBAR", "Sequence"), fn)
  expect_message(s <- topdown:::.readFasta(fn, verbose=TRUE),
                 "Reading sequence from fasta file")
  expect_equal(s, AAString("Sequence"))
  writeLines(c("> FOOBAR", "> Sequence"), fn)
  expect_error(topdown:::.readFasta(fn), "No sequence found")
  unlink(fn)
})

test_that(".readScanHeadsTable", {
  fn <- paste0(tempfile(), ".txt")
  d <- data.frame(MSOrder=c(1, 2, 2, 2, 2),
                  FilterString=c("ms2 100.0001@etd", "ms2 100.0001@hcd",
                                 "ms2 100.0001@hcd", "ms2 100.0007@cid",
                                 "ms2 100.0009@hcd"),
                  Activation1=c("ETD", "ETD", "ETD", "HCD", "CID"),
                  Activation2=c(NA, "HCD", "HCD", "CID", "HCD"),
                  Energy1=c(10, 50, 50, NA, 20),
                  Energy2=c(NA, 30, 30, 20, 10),
                  stringsAsFactors=FALSE)
  write.csv(d, file=fn, row.names=FALSE)
  expect_message(h <- topdown:::.readScanHeadsTable(fn, verbose=TRUE),
                 "Reading 5 header information from file")
  expect_equal(colnames(h),
               c("MSOrder", "FilterString", "Activation1", "Activation2",
                 "Energy1", "Energy2", "Condition",
                 "ETDActivation", "CIDActivation", "HCDActivation",
                 "Activation", "ActivationString", "File"))
  expect_equal(h$MSOrder, rep(2, 4))
  expect_equal(h$ETDActivation, c(50, 50, 0, 0))
  expect_equal(h$CIDActivation, c(0, 0, 20, 20))
  expect_equal(h$HCDActivation, c(30, 30, 0, 10))
  # TODO: FilterStrings are not unique in .experiment.csv files
  # see issue #14
  #expect_equal(h$Condition, c(1, 7, 9))
  expect_equal(h$Condition, c(1, 1:3))
  expect_equal(h$File, rep(gsub("\\.txt$", "", basename(fn)), 4))
  unlink(fn)
})

test_that(".mergeScanConditionAndHeaderInformation", {
  sc <- data.frame(FOO=1:3, Condition=c(1:2, 1), Both=1,
                   File=c("foo", "foo", "bar"))
  hi <- data.frame(BAR=1:5, Condition=c(1, 1, 2, 2, 1), Both=2,
                   File=c("bar", "bar", "bar", "foo", "foo"))
  r <- data.frame(File=c(rep("bar", 3), rep("foo", 2)),
                  Condition=c(1, 1, 2, 1, 2), FOO=c(3, 3, NA, 1, 2),
                  Both.ScanCondition=c(1, 1, NA, 1, 1), BAR=c(1:3, 5:4),
                  Both.HeaderInformation=2)
  expect_equal(topdown:::.mergeScanConditionAndHeaderInformation(sc, hi), r)
})

test_that(".mergeSpectraAndHeaderInformation", {
  fd <- data.frame(x=1:2, Scan=1:2, File="foo", spectrum=1:2)
  hi <- data.frame(File="foo", Scan=1:2, y=3:4, stringsAsFactors=FALSE)
  r <- fd
  r$Scan <- 1:2
  r$y <- 3:4
  r <- r[, c("Scan", "File", "x", "spectrum", "y")]
  expect_equal(topdown:::.mergeSpectraAndHeaderInformation(fd, hi), r)
})
