context("data.frame")

test_that(".aggregateDataFrame", {
    d <- data.frame(id=1:4, col1=1:4, col2=5:8, col3=LETTERS[1:4],
                    stringsAsFactors=FALSE, row.names=letters[1:4])
    r <- data.frame(id=c(1, 3), col1=c(1.5, 3.5), col2=c(5.5, 7.5),
                    col3=LETTERS[c(1, 3)],
                    stringsAsFactors=FALSE, row.names=letters[c(1, 3)])
    expect_equal(topdown:::.aggregateDataFrame(d, rep(1:2, each=2),
                                               ignoreNumCols="id"), r)
    expect_equal(topdown:::.aggregateDataFrame(d, rep(c(10, 2), each=2),
                                               ignoreNumCols="id"), r)
})

test_that(".colsToRle", {
    d <- DataFrame(a=1:10, b=rep(1, 10), c=rep(c("foo", "bar"), each=5))
    r <- DataFrame(a=1:10, b=Rle(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)))
    expect_equal(topdown:::.colsToRle(d), r)
})

test_that(".droplevels", {
    d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(factor(rep(c("foo", "bar"), each=5))),
                   e=factor(rep(1:2, each=5)))
    r <- DataFrame(a=1:5, b=factor(rep(1, 5)),
                   c=Rle(rep("foo", 5)),
                   d=Rle(factor(rep("foo", 5))),
                   e=factor(rep(1, 5)))
    expect_equal(topdown:::.droplevels(d[1:5,]), r)
})

test_that(".dropNonInformativeColumns", {
    d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep("foo", 10)))
    r <- DataFrame(a=1:10,
                   c=Rle(rep(c("foo", "bar"), each=5)))
    expect_equal(topdown:::.dropNonInformativeColumns(d), r)
})

test_that(".orderByColumns", {
    d <- DataFrame(a=10:1, b=c(1:3, 7:1),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep("foo", 10)))
    expect_error(topdown:::.orderByColumns(1:10))
    expect_error(topdown:::.orderByColumns(d, c("foo", "bar")))
    expect_equal(topdown:::.orderByColumns(d, c("c", "a")), 10:1)
})
