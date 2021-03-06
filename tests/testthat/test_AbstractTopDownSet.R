context("AbstractTopDownSet")

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

test_that(".atdsLogMsg", {
    expect_error(.atdsLogMsg(1L, "foo"),
                 "has to be an 'AbstractTopDownSet' object")
    expect_equal(gsub("^\\[[^]]+\\] *", "",
                      .atdsLogMsg(tds, "foobar")@processing),
                 c("Data created.", "foobar; 8 fragments [3;5]."))
    expect_equal(gsub("^\\[[^]]+\\] *", "",
                      .atdsLogMsg(tds, "foobar",
                                            addDim=FALSE)@processing),
                 c("Data created.", "foobar"))
})

test_that("combine", {
    tds$Mz <- 100
    tds$AgcTarget <- 1e5
    tds$EtdReagentTarget <- 1e6
    tds$EtdActivation <- tds$CidActivation <- tds$HcdActivation <- NA_real_
    tds$UvpdActivation <- (1:5) * 1000
    tds$IonInjectionTimeMs <- 1:5
    tds@tolerance <- 5e-6
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

    tds2@rowViews <- FragmentViews("ACE", mass=1:3*90,
                                   type=c("b", "b", "z"),
                                   start=1:3, width=c(1:2, 1),
                                   names=c("b1", "b2", "z1"))
    rownames(tds2@assay) <- c("b1", "b2", "z1")

    tdsr@rowViews <- FragmentViews("ACE", mass=rep(1:3, each=2) * c(90, 100),
                                   type=c("b", "c", "b", "c", "z", "x"),
                                   start=rep(1:3, each=2),
                                   width=rep(c(1:2, 1), each=2),
                                   names=c("b1", "c1", "b2", "c2", "z1", "x1"))
    tdsr@assay <- sparseMatrix(
        i=c(rep(1:4, each=3), rep(5:6, each=2)),
        j=c(2, 4, 6, 1, 3, 5, 2, 4, 10, 1, 3, 9, 4, 8, 3, 7),
        x=c(rep(c(2, 4, 7), 2), rep(c(3, 5, 9), 2), rep(c(6, 8), 2)),
        dimnames=list(c("b1", "c1", "b2", "c2", "z1", "x1"), colnames(a))
    )
    tdsr@processing[7] <- paste("[2017-12-28 15:30:02]",
                                "Combined 8 fragments [3;5] and 8 fragments",
                                "[3;5] into a 16 fragments [6;10]",
                                "TopDownSet object.")

    expect_equal(combine(tds1, tds2), tdsr)
})

test_that("fragmentTypes", {
    expect_error(fragmentType(1L), "doesn't inherit 'AbstractTopDownSet'")
    expect_equal(fragmentType(tds), factor(c("c", "c", "x")))
})

test_that(".inheritsAbstractTopDownSet", {
    expect_true(.inheritsAbstractTopDownSet(new("TopDownSet")))
    expect_true(.inheritsAbstractTopDownSet(new("NCBSet")))
    expect_error(.inheritsAbstractTopDownSet(1L),
                 "doesn't inherit 'AbstractTopDownSet'")
})

test_that(".logdim", {
    expect_error(.logdim(1L))
    expect_equal(.logdim(new("TopDownSet")), "0 fragments [0;0]")
    expect_equal(.logdim(tds), "8 fragments [3;5]")
})

test_that("updateConditionNames", {
    expect_error(updateConditionNames(tds, "Foo"), "must contain names")
    expect_message(updateConditionNames(tds, c("Foo", "Bar", "File"),
                                        verbose=TRUE),
                   "Foo, Bar is/are not present and ignored.")
    expect_message(updateConditionNames(tds, c("File", "Scan"),
                                        verbose=TRUE),
                   "Order of conditions changed.")
    tdn <- tds[, c(4:5, 1:3)]
    colnames(tdn@assay) <- rownames(tdn@colData) <-
        paste0("C", rep(c("bar", "foo"), 2:3), "_", c(4:5, 1:3))
    tdn$Sample <- 1:5
    tdn@processing <- c("[2017-12-23 19:40:00] Data created.",
                        paste0("[2017-12-23 19:40:01] Condition names updated ",
                               "based on: File, Scan. ",
                               "Order of conditions changed. 5 conditions."))
    expect_equal(updateConditionNames(tds, c("File", "Scan"), verbose=FALSE),
                 tdn)
})

test_that("updateMedianInjectionTime", {
    tdn <- tds
    tdn@colData$IonInjectionTimeMs <- 1:5
    tdn@colData$Mz <- tdn@colData$AgcTarget <- rep(1:2, 3:2)
    expect_error(updateMedianInjectionTime(tdn, 1:10), "data length")
    tdr <- tdn
    tdr@colData$MedianIonInjectionTimeMs <- Rle(rep(c(2, 4.5), 3:2))
    tdr@processing <- c(tdr@processing,
                        paste0("[2017-12-20 21:00:00] ",
                               "Recalculate median injection time ",
                               "based on: Mz, AgcTarget."))
    expect_equal(updateMedianInjectionTime(tdn), tdr)
    tdr@colData$Mz <- tdn@colData$Mz <- 1
    tdr@colData$MedianIonInjectionTimeMs <- Rle(3, 5)
    tdr@processing <- c(tdn@processing,
                        paste0("[2017-12-20 21:00:00] ",
                               "Recalculate median injection time."))
    expect_equal(updateMedianInjectionTime(tdn, by=tdn$Mz), tdr)
})
