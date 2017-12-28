context("AbstractTopDownSet")

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

test_that(".atdsLogMsg", {
    expect_error(topdownr:::.atdsLogMsg(1L, "foo"),
                 "has to be an 'AbstractTopDownSet' object")
    expect_equal(gsub("^\\[[^]]+\\] *", "",
                      topdownr:::.atdsLogMsg(tds, "foobar")@processing),
                 c("Data created.", "foobar; 8 fragments [3;5]."))
    expect_equal(gsub("^\\[[^]]+\\] *", "",
                      topdownr:::.atdsLogMsg(tds, "foobar",
                                            addDim=FALSE)@processing),
                 c("Data created.", "foobar"))
})

test_that("fragmentMass", {
    expect_error(fragmentMass(1L), "doesn't inherit 'AbstractTopDownSet'")
    expect_equal(fragmentMass(tds), 1:3 * 100)
})

test_that("fragmentNames", {
    expect_error(fragmentNames(1L), "doesn't inherit 'AbstractTopDownSet'")
    expect_equal(fragmentNames(tds), c("c1", "c2", "x1"))
})

test_that("fragmentTypes", {
    expect_error(fragmentType(1L), "doesn't inherit 'AbstractTopDownSet'")
    expect_equal(fragmentType(tds), factor(c("c", "c", "x")))
})

test_that(".inheritsAbstractTopDownSet", {
    expect_true(topdownr:::.inheritsAbstractTopDownSet(new("TopDownSet")))
    expect_true(topdownr:::.inheritsAbstractTopDownSet(new("NCBSet")))
    expect_error(topdownr:::.inheritsAbstractTopDownSet(1L),
                 "doesn't inherit 'AbstractTopDownSet'")
})

test_that(".logdim", {
    expect_error(topdownr:::.logdim(1L))
    expect_equal(topdownr:::.logdim(new("TopDownSet")), "0 fragments [0;0]")
    expect_equal(topdownr:::.logdim(tds), "8 fragments [3;5]")
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
    expect_equal_TDS(updateConditionNames(tds, c("File", "Scan"),
                                          verbose=FALSE), tdn)
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
    expect_equal_TDS(updateMedianInjectionTime(tdn), tdr)
    tdr@colData$Mz <- tdn@colData$Mz <- 1
    tdr@colData$MedianIonInjectionTimeMs <- Rle(3, 5)
    tdr@processing <- c(tdn@processing,
                        paste0("[2017-12-20 21:00:00] ",
                               "Recalculate median injection time."))
    expect_equal_TDS(updateMedianInjectionTime(tdn, by=tdn$Mz), tdr)
})
