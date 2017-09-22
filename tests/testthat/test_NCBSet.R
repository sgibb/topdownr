context("NCBSet")

expect_equal_NCBSet <- function(object, expected, ..., date=FALSE,
                                info=NULL, label=NULL) {
    if (!date) {
        object@processing <- gsub("^\\[[^]]+\\] *", "", object@processing)
        expected@processing <- gsub("^\\[[^]]+\\] *", "", expected@processing)
    }
    expect_equal(object, expected, ..., info=info, label=label)
}

ncb <- new("NCBSet",
           rowViews=Views(AAString("ACED"), start=1, width=1:3,
                                  names=paste0("bond", 1:3)),
           colData=DataFrame(Scan=1:5, Sample=c(1, 1, 2, 1, 2)),
           assay=sparseMatrix(i=c(1:2, 1:3, 1, 3, 2),
                              j=rep(1:5, c(2, 3, 1, 1, 1)),
                              x=c(1, 2, 1:3, 1, 2, 1)),
           files="aced.fasta",
           processing="[2017-08-19 15:00:00] Data created.")

test_that("[", {

    ncb1 <- new("NCBSet",
                rowViews=Views(AAString("ACED"), start=1, width=1,
                               names="bond1"),
                colData=DataFrame(Scan=1:5, Sample=c(1, 1, 2, 1, 2)),
                assay=sparseMatrix(i=rep(1, 3),
                                   j=1:3,
                                   x=rep(1, 3),
                                   dims=c(1, 5)),
                files="aced.fasta",
                processing=c("[2017-08-19 15:00:00] Data created.",
                             paste0("[2017-08-19 21:48:00] ",
                                    "Subsetted 8 fragments [3;5] to ",
                                    "3 fragments [1;5].")))

    ncb3 <- new("NCBSet",
                rowViews=Views(AAString("ACED"), start=1, width=1:3,
                               names=paste0("bond", 1:3)),
                colData=DataFrame(Scan=1:3, Sample=c(1, 1, 2)),
                assay=sparseMatrix(i=c(1:2, 1:3, 1),
                                   j=rep(1:3, c(2, 3, 1)),
                                   x=c(1, 2, 1:3, 1)),
                files="aced.fasta",
                processing=c("[2017-08-19 15:00:00] Data created.",
                             paste0("[2017-08-19 21:48:00] ",
                                    "Subsetted 8 fragments [3;5] to ",
                                    "6 fragments [3;3].")))

    expect_equal_NCBSet(ncb["bond1"], ncb1)
    expect_equal_NCBSet(ncb[1,], ncb1)
    expect_equal_NCBSet(ncb[,1:3], ncb3)
    expect_warning(ncb[3:1,], "row order")
    expect_warning(ncb[1, drop=TRUE], "'drop' is ignored")
})

test_that("[[ and $", {
    expect_equal(ncb[[1]], 1:5)
    expect_equal(ncb[["Scan"]], 1:5)
    expect_equal(ncb$Scan, 1:5)
    ncbn <- ncb
    ncbn$Scan <- 11:15
    expect_equal(ncbn$Scan, 11:15)
    ncbn[["Scan"]] <- 1:5
    expect_equal(ncbn$Scan, 1:5)
    ncbn$Scan[1:3] <- 11:13
    expect_equal(ncbn$Scan, c(11:13, 4:5))
    ncbn$Scan[1:3] <- 11:13
    expect_equal(ncbn$Scan, c(11:13, 4:5))
})

test_that("accessors", {
    expect_equal(assayData(ncb), ncb@assay)
    expect_equal(colData(ncb), ncb@colData)
    expect_equal(conditionData(ncb), colData(ncb))
    expect_equal(rowViews(ncb), ncb@rowViews)
})

test_that("accessors<-", {
    ncbn <- ncb
    colData(ncbn)[["Scan"]] <- 11:15
    expect_equal(colData(ncbn)$Scan, 11:15)
    ncbn <- ncb
    conditionData(ncbn)[["Scan"]] <- 11:15
    expect_equal(conditionData(ncbn)$Scan, 11:15)
})

test_that("bestConditions", {
    cn <- c("index", "fragments")
    expect_equal(bestConditions(ncb), matrix(c(2, 5, 4, 1), nrow=2,
                                             dimnames=list(NULL, cn)))
    ncbn <- ncb
    ncbn@assay[1, 2] <- 0
    rownames(ncbn@colData) <- LETTERS[1:5]
    expect_equal(bestConditions(ncbn),
                 matrix(c(2:1, 5, 3, 1, 1), nrow=3,
                        dimnames=list(c("B", "A", "E"), cn)))
    expect_equal(bestConditions(ncbn, n=1),
                 matrix(c(2, 3), nrow=1, dimnames=list("B", cn)))
    expect_equal(bestConditions(ncbn, minN=2),
                 matrix(c(2, 3), nrow=1, dimnames=list("B", cn)))
})

test_that("dim", {
    expect_equal(dim(ncb), c(3, 5))
})

test_that("dimnames", {
    expect_equal(dimnames(ncb), list(paste0("bond", 1:3), NULL))
})

test_that("removeEmptyConditions", {
    ncbr <- ncb
    ncbr@assay[cbind(c(1, 3), 3:4)] <- 0L
    ncbr@assay <- drop0(ncbr@assay)
    ncbrr <- ncbr[, c(1:2, 5)]
    ncbrr@processing <- c(ncbrr@processing,
                          paste0("[2017-08-01 22:25:00] ",
                                 "2 empty conditions removed; ",
                                 "6 fragments [3;3]."))
    #ncbr@assay[cbind(c(1, 3), 3:4)] <- 0L
    #ncbr@assay <- drop0(ncbr@assay)
    expect_equal_NCBSet(removeEmptyConditions(ncbr), ncbrr)
})

test_that("show", {
    expect_output(show(new("NCBSet")),
                  "NCBSet object \\([0-9]\\.[0-9]+ Mb\\)")
    expect_output(show(ncb),
                  paste(c("NCBSet object \\([0-9]\\.[0-9]+ Mb\\)",
                          "- - - Protein data - - -",
                          "Amino acid sequence \\(4\\): ACED ",
                          "- - - Fragment data - - -",
                          "Number of N-terminal fragments: 4",
                          "Number of C-terminal fragments: 3",
                          "Number of N- and C-terminal fragments: 1",
                          "- - - Condition data - - -",
                          "Number of conditions: 2 ",
                          "Number of scans: 5 ",
                          "Condition variables \\(2\\): Scan, Sample",
                          "- - - Assay data - - -",
                          "Size of array: 3x5 \\(53\\.33% != 0\\)",
                          "- - - Processing information - - -",
                          "\\[2017-08-19 15:00:00\\] Data created. "),
                        collapse="\n"))
})

test_that("summary", {
    dc <- data.frame(Fragments=c(2, 3, 1, 1, 1), Total=c(3, 6, 1, 2, 1),
                     Min=c(1, 1, 1, 2, 1), Q1=c(1.25, 1.5, 1, 2, 1),
                     Median=c(1.5, 2, 1, 2, 1), Mean=c(1, 2, 1/3, 2/3, 1/3),
                     Q3=c(1.75, 2.5, 1, 2, 1), Max=c(2, 3, 1, 2, 1))
    dr <- data.frame(Fragments=c(3, 3, 2), Total=c(3, 5, 5), Min=c(1, 1, 2),
                     Q1=c(1, 1.5, 2.25), Median=c(1, 2, 2.5), Mean=c(0.6, 1, 1),
                     Q3=c(1, 2, 2.75), Max=1:3)
    expect_equal(summary(ncb), dc)
    expect_equal(summary(ncb, "conditions"), dc)
    expect_equal(summary(ncb, "bonds"), dr)
})

test_that("validObject", {
    ncbn <- ncb
    ncbn@assay <- ncbn@assay[1,,drop=FALSE]
    expect_error(validObject(ncbn), "Mismatch between fragment data")
    ncbn <- ncb
    rownames(ncbn@assay) <- seq_len(nrow(ncbn@assay))
    expect_error(validObject(ncbn), "Mismatch between fragment names")
    ncbn <- ncb
    ncbn@assay <- ncbn@assay[,1,drop=FALSE]
    expect_error(validObject(ncbn), "Mismatch between condition data")
    expect_true(validObject(ncb))
    ncbn <- ncb
    rownames(ncbn@colData) <- seq_len(nrow(ncbn@colData))
    colnames(ncbn@assay) <- rev(seq_len(ncol(ncbn@assay)))
    expect_error(validObject(ncbn), "Mismatch between condition names")
    ncbn <- ncb
    ncbn@assay@x <- ncb@assay@x + 2L
    expect_error(validObject(ncbn), "invalid values")
})
