context("fragments")

test_that(".addAdducts", {
    d <- data.frame(mz=1:6, ion=c("c1", "c2", "c3", "x1", "x2", "y1"),
                    type=rep(c("c", "x", "y"), 3:1), z=1,
                    pos=c(1:3, 1:2, 1),
                    seq=c("A", "AC", "ACE", "E", "CE", "E"),
                    stringsAsFactors=FALSE)
    r <- data.frame(mz=c(2:4, 4), ion=c("cpH1", "cpH2", "cpH3", "ym2H1"),
                    type=c(rep("c", 3), "y"), z=1,
                    pos=c(1:3, 1),
                    seq=c("A", "AC", "ACE", "E"),
                    stringsAsFactors=FALSE)
    a <- data.frame(mass=c(1, -2), name=c("cpH", "ym2H"), to=c("c", "y"))
    expect_error(.addAdducts(d, d),
                 "data.frame must have the columns: 'mass', 'name' and 'to'")
    expect_equal(.addAdducts(d, data.frame()), d)
    expect_equal(.addAdducts(d, a), rbind(d, r))
})

test_that(".addCustomModifications", {
    d <- data.frame(mz=1:6, ion=c("c1", "c2", "c3", "x1", "x2", "y1"),
                    type=rep(c("c", "x", "y"), 3:1), z=1,
                    pos=c(1:3, 1:2, 1),
                    seq=c("A", "AC", "ACE", "E", "CE", "E"),
                    stringsAsFactors=FALSE)
    m <- data.frame(mass=42, name="Acetyl", location="N-term", variable=FALSE)
    r <- d
    r$mz[1:3] <- r$mz[1:3] + 42
    r$type[1:3] <- paste0(r$type[1:3], "_mAcetyl")
    r$ion[1:3] <- paste0(r$ion[1:3], "_mAcetyl")
    expect_error(.addCustomModifications(d, 3, d),
                 "'mass', 'name', 'location' and 'variable'")
    expect_equal(.addCustomModifications(d, 3, data.frame()), d)
    expect_equal(.addCustomModifications(d, 3, m), r)
    m$variable <- TRUE
    expect_equal(.addCustomModifications(d, 3, m), rbind(d, r[1:3,]))
    m <- rbind(m, data.frame(mass=2, name="Foo", location=2, variable=FALSE))
    r <- rbind(d, r[1:3,])
    i <- c(2, 3, 5, 8, 9)
    r$mz[i] <- r$mz[i] + 2
    r$type[i] <- paste0(r$type[i], "_mFoo")
    r$ion[i] <- paste0(r$ion[i], "_mFoo")
    expect_equal(.addCustomModifications(d, 3, m), r)
    r$type[7:9] <- c("c_mAcetyl", "c_mFoo_mAcetyl", "c_mFoo_mAcetyl")
    r$ion[7:9] <- c("c1_mAcetyl", "c2_mFoo_mAcetyl", "c3_mFoo_mAcetyl")
    expect_equal(.addCustomModifications(d, 3, m[2:1,]), r)
})

test_that(".matchFragments", {
    expect_equal(.matchFragments(mz=integer(), fmass=1:3), integer())
    expect_equal(.matchFragments(mz=1:3, fmass=integer()),
                 rep(NA_integer_, 3))
    expect_equal(.matchFragments(c(1, 99, 101), fmass=c(1.1, 100),
                                            redundantFragmentMatch="closest",
                                            tolerance=0.2),
                 as.integer(c(1, 2, NA)))
    expect_equal(.matchFragments(c(1, 98, 101), fmass=c(1.1, 100),
                                            redundantFragmentMatch="closest",
                                            tolerance=0.2),
                 as.integer(c(1, NA, 2)))
    expect_equal(.matchFragments(1:2, 1.5, 1,
                                            redundantFragmentMatch="ignore",
                                            redundantIonMatch="closest",
                                            relative=FALSE),
                 c(1L, 1L))
    expect_equal(.matchFragments(1:2, 1.5, 0.4,
                                            redundantFragmentMatch="closest",
                                            redundantIonMatch="closest",
                                            relative=TRUE),
                 c(1L, NA_integer_))
    expect_equal(.matchFragments(1:2, 1.5, 1,
                                            redundantFragmentMatch="remove",
                                            relative=FALSE),
                 c(NA_integer_, NA_integer_))
    expect_equal(.matchFragments(1:2, 1.5, 1,
                                            redundantFragmentMatch="closest",
                                            relative=FALSE),
                 c(1L, NA_integer_))
    expect_equal(.matchFragments(1.5, 1:2, 1,
                                            redundantIonMatch="remove",
                                            relative=FALSE),
                 NA_integer_)
    expect_equal(.matchFragments(1.5, 1:2, 1,
                                            redundantIonMatch="closest",
                                            relative=FALSE),
                 2L)

    mz <- c(3.4, 3.5, 3.6, 3.7, 8.8, 9.1, 10.9)
    fmass <- 1:10
    tolerance=0.7
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            redundantIonMatch="remove",
                                            redundantFragmentMatch="remove"),
                 .matchFragments(mz, fmass, tolerance))
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            relative=FALSE),
                 rep(NA_integer_, 7))
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            redundantFragmentMatch="closest",
                                            relative=FALSE),
                 as.integer(c(NA, NA, NA, 4, NA, 9, NA)))
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            redundantIonMatch="closest",
                                            relative=FALSE),
                 as.integer(c(3, rep(NA, 6))))
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            redundantFragmentMatch="ignore",
                                            redundantIonMatch="closest",
                                            relative=FALSE),
                 as.integer(c(3, 4, 4, 4, 9, 9, NA)))
    expect_equal(.matchFragments(mz, fmass, tolerance,
                                            redundantFragmentMatch="closest",
                                            redundantIonMatch="closest",
                                            relative=FALSE),
                 as.integer(c(3, NA, NA, 4, NA, 9, NA)))
})

test_that(".reorderSequence", {
    expect_error(.reorderSequence(1:10))
    expect_equal(.reorderSequence("ABCDE"), "ABCDE")
    expect_equal(.reorderSequence("ABCDE", "original"), "ABCDE")
    expect_equal(.reorderSequence("ABCDE", "inverse"), "EDCBA")
    set.seed(2017)
    expect_equal(.reorderSequence("ABCDE", "random"), "DABCE")
})

test_that(".unimod1", {
    d <- data.frame(mz=200:202, pos=1, seq=c("A", "AC", "ACE"),
                    ion=c("c1", "c2", "c3"), stringsAsFactors=FALSE)
    r <- data.frame(mz=200:202 + 42.010565, pos=1, seq=c("A", "AC", "ACE"),
                    ion=c("c1", "c2", "c3"), stringsAsFactors=FALSE)
    expect_equal(.unimod1(d, "ACE"), r)
})

test_that(".unimod4", {
    d <- data.frame(mz=1:5, seq=c("C", "AC", "U", "AE", "AB"),
                    stringsAsFactors=FALSE)
    r <- data.frame(mz=c(57.021464 + 1:3, 4:5), seq=d$seq,
                    stringsAsFactors=FALSE)
    expect_equal(.unimod4(d), r)
})

test_that(".unimod765", {
    expect_equal(.unimod765(c("MACE", "MWE", "EAC")),
                 c("ACE", "MWE", "EAC"))
})

test_that(".cusmod", {
    d <- data.frame(mz=1:6, ion=c("c1", "c2", "c3", "x1", "x2", "y1"),
                    type=rep(c("c", "x", "y"), 3:1), z=1,
                    pos=c(1:3, 1:2, 1),
                    seq=c("A", "AC", "ACE", "E", "CE", "E"),
                    stringsAsFactors=FALSE)
    r <- d
    r$mz[4:6] <- d$mz[4:6] + 0.5
    r$type[4:6] <- paste(d$type[4:6], "mN", sep = "_")
    r$ion[4:6] <- paste(d$ion[4:6], "mN", sep = "_")
    expect_equal(
        .cusmod(d, 3, id = "N", deltamz = 0.5, location = 4, variable = FALSE),
        r
    )

    rvar <- d
    rvar <- rbind(d, r[4:6,], make.row.names = FALSE)
    expect_equal(
        .cusmod(d, 3, id = "N", deltamz = 0.5, location = 4, variable = TRUE),
        rvar
    )
})
