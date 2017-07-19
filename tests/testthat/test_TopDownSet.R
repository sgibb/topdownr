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
           assay=sparseMatrix(i=c(1:2, 1:3, 1, 3, 2),
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
             assay=sparseMatrix(i=rep(1:2, 3),
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
              assay=sparseMatrix(i=rep(1, 3),
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
              assay=sparseMatrix(i=c(1:2, 1:3, 1),
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
  expect_equal_TDS(tds[,Rle(1:3, rep(1, 3))], tdsf)
  expect_warning(tds[3:1,], "row order")
  expect_warning(tds[1, drop=TRUE], "'drop' is ignored")
})

test_that("[[ and $", {
  expect_equal(tds[[1]], 1:5)
  expect_equal(tds[["Scan"]], 1:5)
  expect_equal(tds$Scan, 1:5)
  tdn <- tds
  tdn$Scan <- 11:15
  expect_equal(tdn$Scan, 11:15)
  tdn[["Scan"]] <- 1:5
  expect_equal(tdn$Scan, 1:5)
  tdn$Scan[1:3] <- 11:13
  expect_equal(tdn$Scan, c(11:13, 4:5))
  tdn$Scan[1:3] <- 11:13
  expect_equal(tdn$Scan, c(11:13, 4:5))
})

test_that("accessors", {
  expect_equal(assayData(tds), tds@assay)
  expect_equal(colData(tds), tds@colData)
  expect_equal(conditionData(tds), colData(tds))
  expect_equal(fragmentData(tds), rowViews(tds))
  expect_equal(rowViews(tds), tds@rowViews)
})

test_that("accessors<-", {
  tdn <- tds
  colData(tdn)[["Scan"]] <- 11:15
  expect_equal(colData(tdn)$Scan, 11:15)
  tdn <- tds
  conditionData(tdn)[["Scan"]] <- 11:15
  expect_equal(conditionData(tdn)$Scan, 11:15)
})

test_that("aggregate", {
  tda <- new("TopDownSet",
             rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                    type=c("c", "c", "x"),
                                    start=1:3, width=c(1:2, 1),
                                    names=c("c1", "c2", "x1")),
             colData=DataFrame(Scan=c(1, 4), File=Rle(c("foo", "bar"))),
             assay=sparseMatrix(i=c(1:3, 2:3),
                                 j=rep(1:2, 3:2),
                                 x=c(13/3, 4, 6, 9, 8)),
             files=c("foo.fasta"),
             processing=c("[2017-07-16 14:00:00] Data created.",
                          "[2017-07-16 14:00:01] Aggregated [3;5] to [3;2]."))

  expect_error(aggregate(tds, by="FooBar"), "must be a list")
  expect_equal_TDS(aggregate(tds, by=list(tds$File)), tda)
  expect_equal_TDS(aggregate(tds, by=list(rep(1:2, c(3, 2)))), tda)
})

test_that("dim", {
  expect_equal(dim(tds), c(3, 5))
})

test_that("dimnames", {
  expect_equal(dimnames(tds), list(c("c1", "c2", "x1"), NULL))
})

test_that("filterIntensity", {
  tdl <- tds
  tdl@assay <- sparseMatrix(i=c(1, 3, 2), j=3:5, x=7:9)
  tdl@processing <- c(tdl@processing,
                      "[2017-07-16 14:00:02] 5 intensity values < 7 filtered.")
  expect_error(filterIntensity(tds, 1:2), "length one")
  expect_error(filterIntensity(tds, "c"), "numeric")
  expect_error(filterIntensity(tds, 2), "between 0 and 1")
  expect_error(filterIntensity(tds, -1), "between 0 and 1")
  expect_equal_TDS(filterIntensity(tds, 7, relative=FALSE), tdl)
  tdl@processing[2L] <-
    "[2017-07-16 14:00:02] 5 intensity values < 0.8 (relative) filtered."
  expect_equal_TDS(filterIntensity(tds, 0.8), tdl)
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

test_that(".ncbMap", {
  tds1 <- new("TopDownSet",
              rowViews=FragmentViews("ACE", mass=1:5 * 100,
                                     type=c("c", "c", "x", "y", "z_"),
                                     start=c(1:2, 1, 1, 3),
                                     width=c(1:2, 3, 3, 1),
                                     names=c("c1", "c2", "x3", "y3", "z1_")),
              colData=DataFrame(Scan=1:4),
              assay=sparseMatrix(i=c(1:2, 3:4, 1:4, 5),
                                  j=rep(1:4, c(2, 2, 4, 1)),
                                  x=rep(9, 9)),
              files=c("foo.experiments.csv", "foo.fasta", "foo.mzML",
                      "foo.txt"),
              processing="[2017-07-16 14:00:00] Data created.")

  r <- sparseMatrix(i=rep(1:2, c(3, 4)),
                    j=c(1:3, 1:2, 4:5),
                    x=c(rep(1, 4), 3:1),
                    dims=c(2, 5))

  r1 <- sparseMatrix(i=c(1:2, 1, 1:2),
                     j=c(1, 1, 2, 3, 3),
                     x=c(1, 1, 2, 3, 1),
                     dims=c(3, 4))
  expect_equal(topdown:::.ncbMap(tds), r)
  expect_equal(topdown:::.ncbMap(tds1), r1)
})

test_that("readTopDownFiles", {
  expect_error(readTopDownFiles(".", pattern="FOOBAR"),
               "Could not find any experiments.csv, fasta, mzML, txt files!")
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
                        "Theoretical fragment types \\(2\\): c, x",
                        "Theoretical mass range: \\[100\\.00;300\\.00\\]",
                        "- - - Condition data - - -",
                        "Number of conditions: 5 ",
                        "Condition variables \\(2\\): Scan, File",
                        "- - - Intensity data - - -",
                        "Size of array: 3x5 \\(53\\.33% != 0\\)",
                        "Intensity range: \\[2.00;9.00\\]",
                        "- - - Processing information - - -",
                        "\\[2017-07-16 14:00:00\\] Data created. "),
                      collapse="\n"))
  tdn <- tds
  tdn@assay <- Matrix(0, nrow=3, ncol=5, sparse=TRUE)
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
  tdn@assay <- tdn@assay[1,,drop=FALSE]
  expect_error(validObject(tdn), "Mismatch between fragment data")
  tdn <- tds
  rownames(tdn@assay) <- seq_len(nrow(tdn@assay))
  expect_error(validObject(tdn), "Mismatch between fragment names")
  tdn <- tds
  tdn@assay <- tdn@assay[,1,drop=FALSE]
  expect_error(validObject(tdn), "Mismatch between condition data")
  expect_true(validObject(tds))
  tdn <- tds
  rownames(tdn@colData) <- seq_len(nrow(tdn@colData))
  colnames(tdn@assay) <- rev(seq_len(ncol(tdn@assay)))
  expect_error(validObject(tdn), "Mismatch between condition names")
})
