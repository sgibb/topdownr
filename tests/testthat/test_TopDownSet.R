context("TopDownSet")

expect_equal_TDS <- function(object, expected, ..., date=FALSE,
                            info=NULL, label=NULL) {
  if (!date) {
    object@processing <- gsub("^\\[[^]]+\\] *", "", object@processing)
    expected@processing <- gsub("^\\[[^]]+\\] *", "", expected@processing)
  }
  expect_equal(object, expected, ..., info=info, label=label)
}

tds <- new("TopDownSet",
           rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                  type=c("c", "c", "x"),
                                  start=1:3, width=c(1:2, 1),
                                  names=c("c1", "c2", "x1")),
           colData=DataFrame(Scan=1:5, File=Rle(c(rep("foo", 3),
                                                  "bar", "bar"))),
           assays=sparseMatrix(i=c(1:2, 1:3, 1, 3, 2),
                               j=rep(1:5, c(2, 3, 1, 1, 1)),
                               x=2:9),
           files=c("bar.experiments.csv", "foo.experiments.csv",
                   "foo.fasta", "bar.mzML", "foo.mzML", "bar.txt", "foo.txt"),
           processing="[2017-07-16 14:00:00] Data created.")

test_that("[", {
  tdsc <- new("TopDownSet",
              rowViews=FragmentViews("ACE", mass=1:2 * 100,
                                     type=c("c", "c"),
                                     start=1:2, width=1:2,
                                     names=c("c1", "c2")),
             colData=DataFrame(Scan=1:5, File=Rle(c(rep("foo", 3),
                                                    "bar", "bar"))),
             assays=sparseMatrix(i=rep(1:2, 3),
                                 j=rep(c(1:3, 5), c(2, 2, 1, 1)),
                                 x=c(2:5, 7, 9)),
            files=c("bar.experiments.csv", "foo.experiments.csv",
                    "foo.fasta", "bar.mzML", "foo.mzML", "bar.txt", "foo.txt"),
             processing=c("[2017-07-16 14:00:00] Data created.",
                          "[2017-07-16 14:00:01] Subsetted [3;5] to [2;5]."))
  tdsc1 <- new("TopDownSet",
               rowViews=FragmentViews("ACE", mass=100, type="c",
                                      start=1, width=1, names="c1"),
              colData=DataFrame(Scan=1:5, File=Rle(c(rep("foo", 3),
                                                     "bar", "bar"))),
              assays=sparseMatrix(i=rep(1, 3),
                                  j=1:3,
                                  x=c(2, 4, 7),
                                  dims=c(1, 5)),
             files=c("bar.experiments.csv", "foo.experiments.csv",
                     "foo.fasta", "bar.mzML", "foo.mzML", "bar.txt", "foo.txt"),
              processing=c("[2017-07-16 14:00:00] Data created.",
                           "[2017-07-16 14:00:02] Subsetted [3;5] to [1;5]."))

  tdsf <- new("TopDownSet",
              rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                     type=c("c", "c", "x"),
                                     start=1:3, width=c(1:2, 1),
                                     names=c("c1", "c2", "x1")),
              colData=DataFrame(Scan=1:3, File=Rle(rep("foo", 3))),
              assays=sparseMatrix(i=c(1:2, 1:3, 1),
                                  j=rep(1:3, c(2, 3, 1)),
                                  x=2:7),
              files=c("foo.experiments.csv", "foo.fasta", "foo.mzML",
                      "foo.txt"),
              processing=c("[2017-07-16 14:00:00] Data created.",
                           "[2017-07-16 14:00:03] Subsetted [3;5] to [3;3]."))
  expect_equal_TDS(tds["c"], tdsc)
  expect_equal_TDS(tds["c",], tdsc)
  expect_equal_TDS(tds["c1"], tdsc1)
  expect_equal_TDS(tds[,1:3], tdsf)
  expect_warning(tds[3:1,], "row order")
  expect_warning(tds[1, drop=TRUE], "'drop' is ignored")
})

test_that("dim", {
  expect_equal(dim(tds), c(3, 5))
})

test_that("dimnames", {
  expect_equal(dimnames(tds), list(c("c1", "c2", "x1"), NULL))
})

test_that("fragmentMass", {
  expect_error(fragmentMass(1L), "has to be an 'TopDownSet' object")
  expect_equal(fragmentMass(tds), 1:3 * 100)
})

test_that("fragmentNames", {
  expect_error(fragmentNames(1L), "has to be an 'TopDownSet' object")
  expect_equal(fragmentNames(tds), c("c1", "c2", "x1"))
})

test_that("fragmentTypes", {
  expect_error(fragmentType(1L), "has to be an 'TopDownSet' object")
  expect_equal(fragmentType(tds), factor(c("c", "c", "x")))
})

test_that(".isTopDownSet", {
  expect_true(topdown:::.isTopDownSet(new("TopDownSet")))
  expect_true(topdown:::.isTopDownSet(tds))
  expect_error(topdown:::.isTopDownSet(1L), "has to be an 'TopDownSet' object")
})

test_that("show", {
  expect_output(show(new("TopDownSet")),
                "TopDownSet object \\([0-9]\\.[0-9]+ Mb\\)")
  expect_output(show(tds),
                paste(c("TopDownSet object \\([0-9]\\.[0-9]+ Mb\\)",
                        "- - - Protein data - - -",
                        "Amino acid sequence \\(3\\): ACE ",
                        "- - - Fragment data - - -",
                        "Number of theoretical fragments: 3 ",
                        "Theoretical fragment types \\(2\\): c, x ",
                        "Theoretical mass range: \\[100\\.00;300\\.00\\]",
                        "- - - Condition data - - -",
                        "Number of conditions: 5 ",
                        "Condition variables \\(2\\): Scan, File ",
                        "- - - Intensity data - - -",
                        "Size of array: 3x5 \\(53\\.33% != 0\\)",
                        "Intensity range: \\[2.00;9.00\\]",
                        "- - - Processing information - - -",
                        "\\[2017-07-16 14:00:00\\] Data created. "),
                      collapse="\n"))
  tdn <- tds
  tdn@assays <- Matrix(0, nrow=3, ncol=5, sparse=TRUE)
  expect_output(show(tdn),
                paste(c("Size of array: 3x5 \\(0\\.00% != 0\\)",
                        "Intensity range: \\[NA;NA\\]"),
                      collapse="\n"))
})

test_that(".tdsLogMsg", {
  expect_equal(gsub("^\\[[^]]+\\] *", "",
                    topdown:::.tdsLogMsg(tds, "foobar")@processing),
               c("Data created.", "foobar"))
})

test_that(".validateTopDownSet", {
  tdn <- tds
  tdn@assays <- tdn@assays[1,,drop=FALSE]
  expect_error(validObject(tdn), "Mismatch between fragment data")
  tdn <- tds
  tdn@assays <- tdn@assays[,1,drop=FALSE]
  expect_error(validObject(tdn), "Mismatch between condition data")
  expect_true(validObject(tds))
})
