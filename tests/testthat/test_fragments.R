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
    expect_error(topdown:::.addAdducts(d, d),
                 "data.frame must have the columns: 'mass', 'name' and 'to'")
    expect_equal(topdown:::.addAdducts(d, data.frame()), d)
    expect_equal(topdown:::.addAdducts(d, a), rbind(d, r))
})

test_that(".matchFragments", {
    expect_equal(topdown:::.matchFragments(mz=integer(), fmass=1:3), integer())
    expect_equal(topdown:::.matchFragments(c(1, 99, 101), fmass=c(1.1, 100),
                                           tolerance=0.2),
                 as.integer(c(1, 2, NA)))
    expect_equal(topdown:::.matchFragments(c(1, 98, 101), fmass=c(1.1, 100),
                                           tolerance=0.2),
                 as.integer(c(1, NA, 2)))
})
