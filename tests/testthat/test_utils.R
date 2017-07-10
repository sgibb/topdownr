context("utils")

test_that("cat0", {
  expect_output(topdown:::cat0("foo", "bar"), "foobar")
})

test_that(".filterStringToId", {
  expect_error(topdown:::.filterStringToId(1:3))
  expect_equal(topdown:::.filterStringToId(
    c("FTMS + p NSI sa Full ms2 560.6219@etd50.00@cid7.00 [160.0000-2000.0000]",
      "FTMS + p NSI Full ms2 560.6010@hcd35.00 [160.0000-2000.0000]")),
               c(219L, 10L))
})

test_that(".formatNames", {
  expect_equal(topdown:::.formatNames(c("Monoisotopic M/Z", "SPS Mass 2",
                                        "RT (min)", "MSLevel")),
               c("MonoisotopicMz", "SpsMass2", "RtMin", "MSLevel"))
})

test_that(".fragmentationMethod", {
  d <- expand.grid(ETDActivation=0:1,
                   CIDActivation=0:1,
                   HCDActivation=0:1)
  expect_error(topdown:::.fragmentationMethod(cbind(d, foo=1L)))
  expect_equal(topdown:::.fragmentationMethod(d),
               c("None", "ETD", "CID", "ETcid", "HCD", "EThcd", "HCD/CID", "All"))
})

test_that(".groupBy", {
  x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), na=rep(c(1, NA), 4),
                  stringsAsFactors=FALSE)
  expect_equal(topdown:::.groupBy(x, "LE"), split(x, x$LE))
  expect_equal(topdown:::.groupBy(x, c("ID", "LE")),
               split(x, interaction(as.list(x[, c("ID", "LE")]),
                                    sep=":", lex.order=TRUE)))
  expect_equal(topdown:::.groupBy(x, c("ID", "na")),
               setNames(split(x, x$ID), c("1:1", "2:NA")))
})

test_that(".groupByLabels", {
  x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), na=rep(c(1, NA), 4),
                  stringsAsFactors=FALSE)
  expect_equal(topdown:::.groupByLabels(x, "LE"), x$LE)
  expect_equal(topdown:::.groupByLabels(x, c("ID", "LE")),
               paste(1:2, rep(LETTERS[1:4], each=2), sep=":"))
  expect_equal(topdown:::.groupByLabels(x, c("ID", "na")),
               paste(rep(1:2, 4), rep(c(1, NA), 4), sep=":"))
})

test_that(".massLabel", {
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245)),
               c(750.0001, 1000.8245))
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e5),
               c(750.00001, 1000.80245))
  expect_equal(topdown:::.massLabel(1, 1:999),
               as.double(sprintf("1.%04d", 1:999)))
  expect_error(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e3),
               "at least two digits more than")
})

test_that(".massLabelToId", {
  expect_equal(topdown:::.massLabelToId(c("750.0001", "1000.8245")), c(1, 245))
  expect_equal(topdown:::.massLabelToId(c("750.001", "1000.824"), 2), c(1, 24))
  expect_equal(topdown:::.massLabelToId(sprintf("1000.%04d", 1:999)),
               c(1:999))
})

test_that(".msg", {
  expect_message(topdown:::.msg(TRUE, "foobar"), "foobar")
  expect_message(topdown:::.msg(TRUE, "foo", "bar"), "foobar")
  expect_silent(topdown:::.msg(FALSE, "foobar"))
})

test_that(".nrows", {
  expect_error(topdown:::.nrows(matrix(nrow=2, ncol=2)))
  expect_equal(topdown:::.nrows(list(matrix(nrow=2, ncol=2),
                                     matrix(nrow=3, ncol=2))), 2:3)
})

test_that(".snippet", {
  L <- paste0(LETTERS[1:26], collapse="")
  l <- paste0(letters[1:26], collapse="")
  expect_equal(topdown:::.snippet(L, 100), L)
  expect_equal(topdown:::.snippet(L, 10), "ABCD...XYZ")
  expect_equal(topdown:::.snippet(L, 11), "ABCD...WXYZ")
  expect_equal(topdown:::.snippet(c(l, L), 10), c("abcd...xyz", "ABCD...XYZ"))
  expect_equal(topdown:::.snippet(c(l, L), 11), c("abcd...wxyz", "ABCD...WXYZ"))
})

test_that(".swapFileExt", {
  expect_equal(topdown:::.swapFileExt("foo.xml"), "foo.meth")
  expect_equal(topdown:::.swapFileExt("foo.xml", "bar"), "foo.bar")
})

test_that(".targetedMassListToMz", {
  expect_error(topdown:::.targetedMassListToMz(1:3))
  expect_equal(topdown:::.targetedMassListToMz(c("(mz=1000.12 z=2 name=foo)",
                                                 "(mz=933.99 z=3 name=)")),
               c(1000.1, 933.9))
})

test_that(".topDownFileExtRx", {
  expect_error(topdown:::.topDownFileExtRx("foo"))
  expect_equal(topdown:::.topDownFileExtRx(),
               "\\.experiments\\.csv$|\\.fasta$|\\.mz[Mm][Ll]$|\\.txt$")
  expect_equal(topdown:::.topDownFileExtRx("all"),
               "\\.experiments\\.csv$|\\.fasta$|\\.mz[Mm][Ll]$|\\.raw$|\\.txt$")
  expect_equal(topdown:::.topDownFileExtRx("cfmt"),
               "\\.experiments\\.csv$|\\.fasta$|\\.mz[Mm][Ll]$|\\.txt$")
  expect_equal(topdown:::.topDownFileExtRx("csv"),
               "\\.experiments\\.csv$")
  expect_equal(topdown:::.topDownFileExtRx("mzml"),
               "\\.mz[Mm][Ll]$")
  expect_equal(topdown:::.topDownFileExtRx("txt"),
               "\\.txt$")
})

test_that(".updateAssignmentTableMzId", {
  atab <- data.table(SpectrumId=paste0("F1.S", rep(1:3, each=3)),
                     FragmentId=c(1:3,
                                  3:5,
                                  2:3, 5),
                     MzId=c(1:3,
                            1:3,
                            1:3), key=c("SpectrumId", "FragmentId", "MzId"))
  atab2 <- data.table(SpectrumId=paste0("F1.S", rep(1:3, each=2)),
                      FragmentId=c(1, 3,
                                   3, 5,
                                   2, 5),
                      MzId=c(1:2,
                             1:2,
                             1:2), key=c("SpectrumId", "FragmentId", "MzId"))
  expect_equal(topdown:::.updateAssignmentTableMzId(atab), atab)
  expect_equal(topdown:::.updateAssignmentTableMzId(atab[MzId!=2,]), atab2)
})
