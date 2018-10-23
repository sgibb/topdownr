context("TopDownSet")

tds <- new("TopDownSet",
           rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                  type=c("c", "c", "x"),
                                  start=1:3, width=c(1:2, 1),
                                  names=c("c1", "c2", "x1")),
           colData=DataFrame(Scan=1:5, File=Rle(c(rep("foo", 3),
                                                  "bar", "bar")),
                             row.names=paste("condition",
                                             rep_len(1:3, 5),
                                             rep(1:2, c(3, 2)), sep="_")),
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
                                                      "bar", "bar")),
                                 row.names=paste("condition",
                                                 rep_len(1:3, 5),
                                                 rep(1:2, c(3, 2)), sep="_")),
               assay=sparseMatrix(i=rep(1:2, 3),
                                   j=rep(c(1:3, 5), c(2, 2, 1, 1)),
                                   x=c(2:5, 7, 9)),
              files=c("bar.experiments.csv", "foo.experiments.csv",
                      "foo.fasta", "bar.mzML", "foo.mzML", "bar.txt",
                      "foo.txt"),
               processing=c("[2017-07-16 14:00:00] Data created.",
                            paste0("[2017-07-16 14:00:01] ",
                                   "Subsetted 8 fragments [3;5] to ",
                                   "6 fragments [2;5].")))
    tdsc1 <- new("TopDownSet",
                 rowViews=FragmentViews("ACE", mass=100, type="c",
                                        start=1, width=1, names="c1"),
                colData=DataFrame(Scan=1:5, File=Rle(c(rep("foo", 3),
                                                       "bar", "bar")),
                                  row.names=paste("condition",
                                                  rep_len(1:3, 5),
                                                  rep(1:2, c(3, 2)), sep="_")),
                assay=sparseMatrix(i=rep(1, 3),
                                    j=1:3,
                                    x=c(2, 4, 7),
                                    dims=c(1, 5)),
               files=c("bar.experiments.csv", "foo.experiments.csv",
                       "foo.fasta", "bar.mzML", "foo.mzML", "bar.txt",
                       "foo.txt"),
                processing=c("[2017-07-16 14:00:00] Data created.",
                             paste0("[2017-07-16 14:00:02] ",
                                    "Subsetted 8 fragments [3;5] to ",
                                    "3 fragments [1;5].")))

    tdsf <- new("TopDownSet",
                rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                       type=c("c", "c", "x"),
                                       start=1:3, width=c(1:2, 1),
                                       names=c("c1", "c2", "x1")),
                colData=DataFrame(Scan=1:3, File=Rle(rep("foo", 3)),
                                  row.names=paste("condition", 1:3, rep(1, 3),
                                                  sep="_")),
                assay=sparseMatrix(i=c(1:2, 1:3, 1),
                                    j=rep(1:3, c(2, 3, 1)),
                                    x=2:7),
                files=c("foo.experiments.csv", "foo.fasta", "foo.mzML",
                        "foo.txt"),
                processing=c("[2017-07-16 14:00:00] Data created.",
                             paste0("[2017-07-16 14:00:03] ",
                                    "Subsetted 8 fragments [3;5] to ",
                                    "6 fragments [3;3].")))
    expect_equal(tds["c"], tdsc)
    expect_equal(tds["c",], tdsc)
    expect_equal(tds["c1"], tdsc1)
    expect_equal(tds[,1:3], tdsf)
    expect_equal(tds[,Rle(1:3, rep(1, 3))], tdsf)
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
    tds$Sample <- c(2, 2, 2, 10, 10)
    tda <- new("TopDownSet",
               rowViews=FragmentViews("ACE", mass=1:3 * 100,
                                      type=c("c", "c", "x"),
                                      start=1:3, width=c(1:2, 1),
                                      names=c("c1", "c2", "x1")),
               colData=DataFrame(Scan=c(1, 4), File=Rle(c("foo", "bar")),
                                 Sample=c(2, 10),
                                 row.names=paste0("condition_1_", 1:2)),
               assay=sparseMatrix(i=c(1:3, 2:3),
                                   j=rep(1:2, 3:2),
                                   x=c(13/3, 4, 6, 9, 8)),
               files=tds@files,
               processing=c("[2017-07-16 14:00:00] Data created.",
                            paste0("[2017-07-16 14:00:01] Aggregated ",
                                   "8 fragments [3;5] to ",
                                   "5 fragments [3;2].")))
    expect_error(aggregate(tds, by="FooBar"), "same length")
    expect_equal(aggregate(tds, by=list(tds$File)), tda)
    expect_equal(aggregate(tds, by=list(rep(1:2, c(3, 2)))), tda)
    expect_equal(aggregate(tds, by=list(rep(c(2, 10), c(3, 2)))), tda)
})

test_that("combine", {
    tds$Mz <- 100
    tds$AgcTarget <- 1e5
    tds$EtdReagentTarget <- 1e6
    tds$EtdActivation <- tds$CidActivation <- tds$HcdActivation <- NA_real_
    tds$UvpdActivation <- (1:5) * 1000
    tds$IonInjectionTimeMs <- 1:5
    tds@tolerance <- 5e-6
    tds@redundantMatching <- c("remove", "remove")
    rownames(tds@assay) <- rownames(tds)
    tds <- updateConditionNames(tds)
    tds1 <- tds2 <- tds
    o <- c(matrix(1:10, nrow=2, byrow=TRUE))
    a <- cbind(tds1@assay, tds2@assay)[, o]
    cd <- .colsToRle(rbind(tds1@colData, tds2@colData)[o, ])
    colnames(a) <- rownames(cd) <- paste0("C", rep((1:5) * 1000, each=2), "_",
                                          rep(1:2, 5))
    tdsr <- new("TopDownSet",
                rowViews=tds@rowViews,
                colData=cd,
                assay=a,
                tolerance=5e-6,
                files=unique(tds1@files, tds2@files),
                processing=c(tds1@processing, tds2@processing,
                             paste("[2017-12-28 15:30:00]",
                                   "Condition names updated based on: Mz,",
                                   "AgcTarget, EtdReagentTarget,",
                                   "EtdActivation, CidActivation,",
                                   "HcdActivation, UvpdActivation. Order of",
                                   "conditions changed. 5 conditions."),
                             paste("[2017-12-28 15:30:01]",
                                   "Recalculate median injection time based",
                                   "on: Mz, AgcTarget."),
                             paste("[2017-12-28 15:30:02]",
                                   "Combined 8 fragments [3;5] and 8 fragments",
                                   "[3;5] into a 16 fragments [3;10]",
                                   "TopDownSet object.")))
    tdsr$MedianIonInjectionTimeMs <- Rle(3, 10)
    expect_equal(combine(tds1, tds2), tdsr)
    tds2@tolerance <- 1e-6
    expect_warning(combine(tds1, tds2), "Matching tolerance differs")
    expect_equal(suppressWarnings(combine(tds1, tds2)), tdsr)
    tds2@tolerance <- tds1@tolerance
    tds2@redundantMatching <- c("closest", "closest")
    expect_warning(combine(tds1, tds2), "Matching strategies differ")
    expect_equal(suppressWarnings(combine(tds1, tds2)), tdsr)
})

test_that("conditionNames", {
    expect_equal(conditionNames(tds),
                 paste("condition",
                       rep_len(1:3, 5),
                       rep(1:2, c(3, 2)), sep="_"))
})

test_that("condition2data.frame", {
    expect_error(topdownr:::.condition2data.frame(1:10),
                 "has to be an 'TopDownSet' object")
    expect_error(topdownr:::.condition2data.frame(tds),
                 "more than one condition")
    skip_if_not_installed("topdownrdata", "0.2")
    suppressWarnings(tds <- readTopDownFiles(
        topdownrdata::topDownDataPath("myoglobin"),
        pattern=".*fasta.gz$|1211_.*1e6_1",
        neutralLoss=NULL,
        tolerance=25e-6
    ))
    s <- topdownr:::.readSpectrum(tds@files[grep("mzML.gz", tds@files)],
                                  tds$SpectrumIndex[3])
    d <- data.frame(
        mz=s[,1], intensity=s[,2], fragment="",
        type=factor(
            "None",
            levels=c("None", "N-terminal", "C-terminal", "Bidirectional"),
            ordered=TRUE
        ),
        stringsAsFactors=FALSE
    )
    i <- c(29, 34, 42, 45, 60, 76, 77, 78, 84, 85, 95, 97, 98, 99, 100, 114, 117)
    d$fragment[i] <- c("y48", "y51", "y58", "x60", "y75", "z90", "b89", "c89",
        "y94", "x94", "y101", "y105", "z106", "y106", "y112", "y152", "c152")
    d$type[i] <- paste0(c("C", "C", "C", "C", "C", "C", "N", "N", "C", "C",
        "C", "C", "C", "C", "C", "C", "N"), "-terminal")
    expect_equal(topdownr:::.condition2data.frame(tds[, 3]), d)
    s <- topdownr:::.readSpectrum(tds@files[grep("mzML.gz", tds@files)],
                                  tds$SpectrumIndex[103])
    d <- data.frame(
        mz=s[,1], intensity=s[,2], fragment="",
        type=factor(
            "None",
            levels=c("None", "N-terminal", "C-terminal", "Bidirectional"),
            ordered=TRUE
        ),
        stringsAsFactors=FALSE
    )
    i <- c(23, 26)
    d$fragment[i] <- c("z152", "b152")
    d$type[i] <- paste0(c("C", "N"), "-terminal")
    expect_equal(topdownr:::.condition2data.frame(tds[, 103]), d)
})

test_that("dim", {
    expect_equal(dim(tds), c(3, 5))
})

test_that("dimnames", {
    expect_equal(dimnames(tds), list(c("c1", "c2", "x1"),
                                     paste("condition",
                                           rep_len(1:3, 5),
                                           rep(1:2, c(3, 2)), sep="_")))
})

test_that("filterCv", {
    tdfit <- tds
    tdfit$Sample <- rep(1:2, 3:2)

    expect_error(filterCv(tdfit, threshold="A"),
                 "has to be a 'numeric'")
    expect_error(filterCv(tdfit, threshold=1:2),
                 "length one")
    expect_error(filterCv(tdfit, threshold=-1),
                 "greater than 0")

    expect_equal(filterCv(tdfit, threshold=100), tdfit)
    tdfitr <- tdfit
    tdfitr@assay[c(1, 4, 7)] <- 0L
    tdfitr@assay <- drop0(tdfitr@assay)
    tdfitr@processing <- c(tdfitr@processing,
                           paste0("[2017-08-04 18:05:00] 3 fragments with ",
                                  "CV > 40% filtered; 5 fragments [3;5]."))
    expect_equal(filterCv(tdfit, threshold=40), tdfitr)
})

test_that("filterInjectionTime", {
    tdfit <- tds
    tdfit$IonInjectionTimeMs <- c(1:3, 2:3)
    tdfit$MedianIonInjectionTimeMs <- rep(1:2, 3:2)
    tdfit$Sample <- rep(1:2, 3:2)

    expect_error(filterInjectionTime(tdfit, maxDeviation="A"),
                 "has to be a 'numeric'")
    expect_error(filterInjectionTime(tdfit, maxDeviation=1:2),
                 "length one")
    expect_error(filterInjectionTime(tdfit, maxDeviation=-1),
                 "greater than 0")

    expect_equal(filterInjectionTime(tdfit, maxDeviation=100, keepTopN=5),
                 tdfit)
    tdfitr <- tdfit[, c(1, 4)]
    tdfitr@processing <- c(tdfitr@processing,
                           paste0("[2017-07-28 16:00:02] 3 scans filtered ",
                                  "with injection time deviation >= 0.5 or ",
                                  "rank >= 6; 3 fragments [3;2]."))
    expect_equal(filterInjectionTime(tdfit, maxDeviation=0.5, keepTopN=5),
                 tdfitr)
    tdfitr <- tdfit[, -3]
    tdfitr@processing <- c(tdfitr@processing,
                           paste0("[2017-07-28 16:00:02] 1 scan filtered ",
                                  "with injection time deviation >= 5 or ",
                                  "rank >= 3; 7 fragments [3;4]."))
    expect_equal(filterInjectionTime(tdfit, maxDeviation=5, keepTopN=2),
                 tdfitr)
    tdfitr <- tdfit[, -(2:3)]
    tdfitr@processing <- c(tdfitr@processing,
                           paste0("[2017-07-28 16:00:02] 2 scans filtered ",
                                  "with injection time deviation >= 0.6 or ",
                                  "rank >= 3; 4 fragments [3;3]."))
    expect_equal(filterInjectionTime(tdfit, maxDeviation=0.6, keepTopN=2),
                 tdfitr)
})

test_that("filterIntensity", {
   tdl <- tds
   tdl@assay <- sparseMatrix(i=c(1, 3, 2), j=3:5, x=7:9)
   tdl@processing <- c(tdl@processing,
                       paste0("[2017-07-16 14:00:02] ",
                              "5 intensity values < 7 filtered; ",
                              "3 fragments [3;5]."))
   expect_error(filterIntensity(tds, 1:2), "length one")
   expect_error(filterIntensity(tds, "c"), "numeric")
   expect_error(filterIntensity(tds, 2), "between 0 and 1")
   expect_error(filterIntensity(tds, -1), "between 0 and 1")
   expect_equal(filterIntensity(tds, 7, relative=FALSE), tdl)
   tdl@processing[2L] <- paste0("[2017-07-16 14:00:02] ",
                                "5 intensity values < 0.8 (relative) filtered; ",
                                "3 fragments [3;5].")
   expect_equal(filterIntensity(tds, 0.8), tdl)
})

test_that("filterNonReplicatedFragments", {
    tdfit <- tds
    tdfit$Sample <- rep(1:2, 3:2)

    expect_error(filterNonReplicatedFragments(tdfit, minN="A"),
                 "has to be a 'numeric'")
    expect_error(filterNonReplicatedFragments(tdfit, minN=1:2),
                 "length one")

    expect_equal(filterNonReplicatedFragments(tdfit, minN=0),
                 tdfit)
    tdfitr <- tdfit
    tdfitr@assay <- sparseMatrix(i=c(1:2, 1:2, 1),
                                 j=rep(1:3, c(2, 2, 1)),
                                 x=c(2:5, 7), dims=c(3, 5))
    tdfitr@processing <- c(tdfitr@processing,
                           paste0("[2017-07-31 22:15:00] 3 intensity values ",
                                  "of fragments replicated < 2 times ",
                                  "filtered; 5 fragments [3;5]."))
    expect_equal(filterNonReplicatedFragments(tdfit, minN=2), tdfitr)
})

test_that(".isTopDownSet", {
    expect_true(topdownr:::.isTopDownSet(new("TopDownSet")))
    expect_true(topdownr:::.isTopDownSet(tds))
    expect_error(topdownr:::.isTopDownSet(1L),
                 "has to be an 'TopDownSet' object")
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
                      dims=c(2, 5),
                      dimnames=list(paste0("bond", 1:2), colnames(tds)))

    r1 <- sparseMatrix(i=c(1:2, 1, 1:2),
                       j=c(1, 1, 2, 3, 3),
                       x=c(1, 1, 2, 3, 1),
                       dims=c(3, 4),
                       dimnames=list(paste0("bond", 1:3), c()))

    expect_equal(topdownr:::.ncbMap(tds), r)
    expect_equal(topdownr:::.ncbMap(tds1), r1)
})

test_that("normalize", {
    tds$TotIonCurrent <- c(10, 20, 10, 10, 20)
    tdn <- tds
    tdn@assay <- t(t(tds@assay) / tdn$TotIonCurrent)
    tdn@processing <- c(tdn@processing,
                        paste0("[2017-08-06 14:50:00] ",
                               "Intensity values normalized ",
                               "to TIC."))
    expect_equal(normalize(tds, method="TIC"), tdn)
})

test_that("readTopDownFiles", {
    expect_error(readTopDownFiles(".", pattern="FOOBAR"),
                 "Could not find any csv, fasta, mzML, txt files!")
})

test_that("removeEmptyConditions", {
    tdr <- tds
    tdr@assay[cbind(c(1, 3), 3:4)] <- 0L
    tdr@assay <- drop0(tdr@assay)
    tdrr <- tdr[, c(1:2, 5)]
    tdrr@processing <- c(tdrr@processing,
                         paste0("[2017-08-01 22:25:00] ",
                                "2 empty conditions removed; ",
                                "6 fragments [3;3]."))
    expect_equal(removeEmptyConditions(tdr), tdrr)
})

test_that("show", {
    expect_output(show(new("TopDownSet")),
                  "TopDownSet object \\([0-9]\\.[0-9]+ Mb\\)")
    tds$Sample <- rep(1:2, 3:2)
    expect_output(show(tds),
                  paste(c("TopDownSet object \\([0-9]\\.[0-9]+ Mb\\)",
                          "- - - Protein data - - -",
                          "Amino acid sequence \\(3\\): ACE ",
                          "- - - Fragment data - - -",
                          "Number of theoretical fragments: 3 ",
                          "Theoretical fragment types \\(2\\): c, x",
                          "Theoretical mass range: \\[100\\.00;300\\.00\\]",
                          "- - - Condition data - - -",
                          "Number of conditions: 2 ",
                          "Number of scans: 5 ",
                          "Condition variables \\(3\\): Scan, File, Sample",
                          "- - - Intensity data - - -",
                          "Size of array: 3x5 \\(53\\.33% != 0\\)",
                          "Number of matched fragments: 8 ",
                          "Intensity range: \\[2.00;9.00\\]",
                          "- - - Processing information - - -",
                          "\\[2017-07-16 14:00:00\\] Data created. "),
                        collapse="\n"))
    tdn <- tds
    tdn@assay <- Matrix(0, nrow=3, ncol=5, sparse=TRUE)
    expect_output(show(tdn),
                  paste(c("Size of array: 3x5 \\(0\\.00% != 0\\)",
                          "Number of matched fragments: 0 ",
                          "Intensity range: \\[NA;NA\\]"),
                        collapse="\n"))
    tdn@rowViews@metadata$modifications <-
        c("Carbamidomethyl", "Met-loss+Acetyl")
    expect_output(show(tdn),
                  "Modifications \\(2\\): Carbamidomethyl, Met-loss\\+Acetyl")
})

test_that("summary", {
    dc <- data.frame(Fragments=c(2, 3, 1, 1, 1), Total=c(5, 15, 7:9),
                     Min=c(2, 4, 7:9), Q1=c(2.25, 4.5, 7:9),
                     Median=c(2.5, 5.0, 7:9), Mean=c(5/3, 5, 7/3, 8/3, 3),
                     Q3=c(2.75, 5.5, 7:9), Max=c(3, 6, 7:9))
    dr <- data.frame(Fragments=c(3, 3, 2), Total=c(13, 17, 14), Min=c(2:3, 6),
                     Q1=c(3:4, 6.5), Median=c(4:5, 7), Mean=c(2.6, 3.4, 2.8),
                     Q3=c(5.5, 7.0, 7.5), Max=c(7, 9, 8))
    expect_equal(summary(tds), dc)
    expect_equal(summary(tds, "conditions"), dc)
    expect_equal(summary(tds, "fragments"), dr)
})

test_that("validObject", {
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

test_that("as(\"MSnSet\")", {
    msn <- new("MSnSet",
               exprs=matrix(c(2:3, 0, 4:6, 7, rep(0, 4), 8, 0, 9, 0), nrow=3,
                            dimnames=dimnames(tds)),
               phenoData=as(data.frame(Scan=1:5,
                                       File=rep(c("foo", "bar"), 3:2),
                                       row.names=colnames(tds),
                                       stringsAsFactors=FALSE),
                            "AnnotatedDataFrame"),
               featureData=as(data.frame(fragment=c("A", "CE", "E"),
                                         start=1:3, end=c(1, 3, 3), width=c(1, 2, 1),
                                         name=c("c1", "c2", "x1"),
                                         type=factor(c("c", "c", "x")),
                                         mass=(1:3)*100, z=1,
                                         row.names=c("c1", "c2", "x1")),
                              "AnnotatedDataFrame"))
    msn@processingData <- new("MSnProcess", files=tds@files,
                              processing="[2017-07-16 14:00:00] Data created.")
    expect_equal(as(tds, "MSnSet"), msn)
})

test_that("as(\"NCBSet\")", {
    ncb <- new("NCBSet",
               rowViews=Views(AAString("ACE"), start=1, width=1:2,
                              names=paste0("bond", 1:2)),
               colData=DataFrame(tds@colData, AssignedIntensity=c(5, 15, 7:9),
                                 row.names=paste("condition",
                                                 rep_len(1:3, 5),
                                                 rep(1:2, c(3, 2)), sep="_")),
               assay=sparseMatrix(i=rep(1:2, c(3, 4)),
                                  j=c(1:3, 1:2, 4:5),
                                  x=c(rep(1, 4), 3, 2, 1),
                                  dimnames=list(c("bond1", "bond2"),
                                                paste("condition",
                                                      rep_len(1:3, 5),
                                                      rep(1:2, c(3, 2)),
                                                      sep="_"))),
               files=tds@files,
               processing=c(tds@processing,
                            paste("[2017-08-20 16:30:00]",
                                  "Coerced TopDownSet into an NCBSet object;",
                                  "7 fragments [2;5].")))
    expect_equal(as(tds, "NCBSet"), ncb)
})
