context("FragmentViews")

fv <- FragmentViews("ACE", start=3:1, end=3, names=paste0("x", 1:3),
                    mass=c(174.039696, 334.070346, 405.107456),
                    type=rep("x", 3), z=rep(1, 3),
                    metadata=list(modifications=c("Carbamidomethyl",
                                                  "Met-loss+Acetyl")))

test_that("constructor", {
    expect_error(topdown:::.calculateFragments("ACE", modifications="FOO"))
    expect_equal(topdown:::.calculateFragments("ACE", type="x",
                                               neutralLoss=list(water=NULL)),
                 fv)
    fv1 <- FragmentViews("AC", start=1:2, end=2, names=paste0("x", 2:1),
                         mass=c(276.064870 + 42.010565, 205.027760),
                         type=rep("x", 2), z=rep(1, 2),
                         metadata=list(modifications=c("Carbamidomethyl",
                                                       "Met-loss+Acetyl")))
    expect_equal(topdown:::.calculateFragments("MAC", type="x",
                                               neutralLoss=list(water=NULL)),
                 fv1, tolerance=1e-7)
    fv2 <- FragmentViews("AC", start=1:2, end=2, names=paste0("x", 2:1),
                         mass=c(276.064870 + 42.010565, 205.027760) - 57.021464,
                         type=rep("x", 2), z=rep(1, 2),
                         metadata=list(modifications=c("Met-loss+Acetyl")))
    expect_equal(topdown:::.calculateFragments("MAC", type="x",
                                               modifications="Met-loss+Acetyl",
                                               neutralLoss=list(water=NULL)),
                 fv2, tolerance=1e-7)
    fv3 <- FragmentViews("CA", start=2:1, end=2, names=paste0("x", 1:2),
                         mass=c(116.034216, 276.064866 - 57.021464),
                         type=rep("x", 2), z=rep(1, 2),
                         metadata=list(modifications=c("Met-loss+Acetyl")))
    expect_equal(topdown:::.calculateFragments("MCA", type="x",
                                               modifications="Met-loss+Acetyl",
                                               neutralLoss=list(water=NULL)),
                 fv3, tolerance=1e-7)
})

test_that("show", {
    expect_output(show(fv),
                  paste(c("FragmentViews on a 3-letter sequence:",
                          "  ACE",
                          "Modifications:",
                          "  Carbamidomethyl",
                          "  Met-loss\\+Acetyl",
                          "Views:",
                          "    start end width   mass type z *",
                          "\\[1\\]     3   3     1 174\\.04 x    1 \\[E\\] *",
                          "\\[2\\]     2   3     2 334\\.07 x    1 \\[CE\\] *",
                          "\\[3\\]     1   3     3 405\\.11 x    1 \\[ACE\\]"),
                        collapse="\n"))
})

test_that("validity", {
    expect_true(validObject(fv))
    d <- elementMetadata(fv)
    elementMetadata(fv) <- NULL
    expect_error(validObject(fv), "'mass', 'type', 'z' are missing")
    elementMetadata(fv) <- d
    elementMetadata(fv)$z <- NULL
    expect_error(validObject(fv), "'z' is missing")
    elementMetadata(fv) <- d
    elementMetadata(fv)$mass <- 1L:3L
    expect_error(validObject(fv), "'mass' has to be of type double")
    elementMetadata(fv) <- d
    elementMetadata(fv)$mass <- as.double(3:1)
    expect_error(validObject(fv), "'mass' has to be sorted")
})
