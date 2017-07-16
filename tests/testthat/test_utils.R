context("utils")

test_that("cat0", {
  expect_output(topdown:::cat0("foo", "bar"), "foobar")
})

test_that(".colsToRle", {
  d <- DataFrame(a=1:10, b=rep(1, 10), c=rep(c("foo", "bar"), each=5))
  r <- DataFrame(a=1:10, b=Rle(rep(1, 10)), c=Rle(rep(c("foo", "bar"), each=5)))
  expect_equal(topdown:::.colsToRle(d), r)
})

test_that(".droplevels", {
  d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                 c=Rle(rep(c("foo", "bar"), each=5)),
                 d=Rle(factor(rep(c("foo", "bar"), each=5))),
                 e=factor(rep(1:2, each=5)))
  r <- DataFrame(a=1:5, b=factor(rep(1, 5)),
                 c=Rle(rep("foo", 5)),
                 d=Rle(factor(rep("foo", 5))),
                 e=factor(rep(1, 5)))
  expect_equal(topdown:::.droplevels(d[1:5,]), r)
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
                                        "RT (min)", "MSLevel",
                                        "Multi.Inject.Info", "RF.Comp...ppm")),
               c("MonoisotopicMZ", "SpsMass2", "RtMin", "MSLevel",
                 "MultiInjectInfo", "RfCompPpm"))
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

test_that(".logmsg", {
  expect_equal(topdown:::.logmsg("foo"),
               paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] foo"))
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

test_that(".subset", {
  expect_error(topdown:::.subset(1:2, 10, letters[1:2]))
  expect_error(topdown:::.subset(c(1, NA, 2), 10, letters[1:10]),
               "'NA' is not supported")
  expect_error(topdown:::.subset(list(foo=1:10), 10, letters[1:10]),
               "Unknown")
  expect_equal(topdown:::.subset(1:2, 10, letters[1:10]), 1:2)
  expect_equal(topdown:::.subset(c(TRUE, TRUE, rep(FALSE, 8)), 10,
                                 letters[1:10]), 1:2)
  expect_equal(topdown:::.subset(c("a", "b"), 10, letters[1:10]), 1:2)
})

test_that(".subsetByCharacter", {
  expect_error(topdown:::.subsetByCharacter(1:2, LETTERS[1:2]))
  expect_error(topdown:::.subsetByCharacter(letters[1:2], TRUE))
  expect_error(topdown:::.subsetByCharacter(letters[1:2], LETTERS[1:2]),
               "Subscript out of bound: 'a', 'b'")
  expect_equal(topdown:::.subsetByCharacter(letters[1:2], letters[4:1]), 4:3)
  expect_equal(topdown:::.subsetByCharacter(letters[1:2]), integer())
})

test_that(".subsetByLogical", {
  expect_error(topdown:::.subsetByLogical(1:2, 10))
  expect_error(topdown:::.subsetByLogical(TRUE, TRUE))
  expect_error(topdown:::.subsetByLogical("foo", 10))
  expect_equal(topdown:::.subsetByLogical(TRUE, 10), 1:10)
  expect_equal(topdown:::.subsetByLogical(c(TRUE, FALSE), 10),
               seq(1, 10, by=2))
  expect_equal(topdown:::.subsetByLogical(rep(TRUE, 10), 10), 1:10)
  expect_equal(topdown:::.subsetByLogical(rep(TRUE, 12), 10), 1:10)
})

test_that(".subsetByNumeric", {
  expect_error(topdown:::.subsetByNumeric(TRUE, 10))
  expect_error(topdown:::.subsetByNumeric(1:10, TRUE))
  expect_error(topdown:::.subsetByNumeric("foo", 10))
  expect_error(topdown:::.subsetByNumeric(c(-1, 3, 12), 10),
               "Subscript out of bound: '-1', '12'")
  expect_equal(topdown:::.subsetByNumeric(1:10, 20), 1:10)
})

test_that(".subsetFiles", {
  expect_equal(topdown:::.subsetFiles(
                c("foo.experiments.csv", "foo.mzML", "bar.txt"), "foo"),
               c(TRUE, TRUE, FALSE))
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
