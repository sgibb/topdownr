context("FragmentViews")

fv <- FragmentViews("AACE", start=1, end=1:3, names=paste0("b", 1:3),
                    mass=c(114.054951, 185.092061, 345.122711),
                    type=rep("b", 3), z=rep(1, 3),
                    metadata=list(modifications=c("Carbamidomethyl",
                                                  "Acetyl", "Met-loss"),
                                  mass=473.158029))

test_that("constructor", {
    expect_error(topdown:::.calculateFragments("AACE", modifications="FOO"))
    expect_equal(topdown:::.calculateFragments("MAACE", type="b"), fv)

    ## default: Acetylation, Carboxyamidomethylation and Met-loss
    expect_equal(topdown:::.calculateFragments("MAACE", type="b"),
                 fv)

    ## Acetylation but no Carboxyamidomethylation/Met-loss
    fv2 <- fv
    fv2@metadata$mass <- fv@metadata$mass - 57.021464
    fv2@elementMetadata$mass[3] <- fv@elementMetadata$mass[3] - 57.021464
    fv2@metadata$modifications <- c("Acetyl")
    expect_equal(topdown:::.calculateFragments("AACE", type="b",
                                               modifications="Acetyl"),
                 fv2)
    ## just Met-loss without Acetylation
    fv3 <- fv
    fv3@metadata$mass <- fv@metadata$mass - 57.021464 - 42.010565
    fv3@elementMetadata$mass <- fv@elementMetadata$mass - 42.010565
    fv3@elementMetadata$mass[3] <- fv3@elementMetadata$mass[3] - 57.021464
    fv3@metadata$modifications <- c("Met-loss")
    expect_equal(topdown:::.calculateFragments("MAACE", type="b",
                                               modifications="Met-loss"),
                 fv3)
    ## no modification
    fv4 <- fv
    fv4@elementMetadata$mass <- fv@elementMetadata$mass - 42.010565
    fv4@elementMetadata$mass[3] <- fv4@elementMetadata$mass[3] - 57.021464
    fv4@metadata <- list(modifications=NULL,
                         mass=fv@metadata$mass - 57.021464 - 42.010565)
    expect_equal(topdown:::.calculateFragments("AACE", type="b",
                                               modifications=NULL), fv4)
})

test_that("show", {
    expect_output(show(fv),
                  paste(c("FragmentViews on a 4-letter sequence:",
                          "  AACE",
                          "Mass:",
                          "  473.158029",
                          "Modifications:",
                          "  Carbamidomethyl",
                          "  Acetyl",
                          "  Met-loss",
                          "Views:",
                          "    start end width   mass name type z *",
                          "\\[1\\]     1   1     1 114\\.05 b1   b    1 \\[A\\] *",
                          "\\[2\\]     1   2     2 185\\.09 b2   b    1 \\[AA\\] *",
                          "\\[3\\]     1   3     3 345\\.12 b3   b    1 \\[AAC\\]"),
                        collapse="\n"))
    fv@metadata <- list()
    expect_output(show(fv),
                  paste(c("FragmentViews on a 4-letter sequence:",
                          "  AACE",
                          "Views:",
                          "    start end width   mass name type z *",
                          "\\[1\\]     1   1     1 114\\.05 b1   b    1 \\[A\\] *",
                          "\\[2\\]     1   2     2 185\\.09 b2   b    1 \\[AA\\] *",
                          "\\[3\\]     1   3     3 345\\.12 b3   b    1 \\[AAC\\]"),
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

test_that("as(\"data.frame\")", {
    expect_equal(as(fv, "data.frame"),
                 data.frame(fragment=c("A", "AA", "AAC"),
                            start=1,
                            end=1:3,
                            width=1:3,
                            name=paste0("b", 1:3),
                            type=factor(rep("b", 3)),
                            mass=c(114.054951, 185.092061, 345.122711),
                            z=1,
                            row.names=paste0("b", 1:3)))
})
