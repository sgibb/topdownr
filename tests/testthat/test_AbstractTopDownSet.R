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
    expect_equal(gsub("^\\[[^]]+\\] *", "",
                      topdown:::.atdsLogMsg(tds, "foobar")@processing),
                 c("Data created.", "foobar"))
})

test_that(".inheritsAbstractTopDownSet", {
    expect_true(topdown:::.inheritsAbstractTopDownSet(new("TopDownSet")))
    expect_true(topdown:::.inheritsAbstractTopDownSet(new("NCBSet")))
    expect_error(topdown:::.inheritsAbstractTopDownSet(1L),
                 "doesn't inherit 'AbstractTopDownSet'")
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


