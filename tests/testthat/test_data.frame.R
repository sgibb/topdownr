context("data.frame")

test_that(".aggregateDataFrame", {
    d <- data.frame(id=1:4, col1=1:4, col2=5:8, col3=LETTERS[1:4],
                    stringsAsFactors=FALSE, row.names=letters[1:4])
    r <- data.frame(id=c(1, 3), col1=c(1.5, 3.5), col2=c(5.5, 7.5),
                    col3=LETTERS[c(1, 3)],
                    stringsAsFactors=FALSE, row.names=letters[c(1, 3)])
    expect_equal(topdownr:::.aggregateDataFrame(d, rep(1:2, each=2),
                                               ignoreNumCols="id"), r)
    expect_equal(topdownr:::.aggregateDataFrame(d, rep(c(10, 2), each=2),
                                               ignoreNumCols="id"), r)
})

test_that(".colsToRle", {
    d <- DataFrame(a=1:10, b=rep(1, 10), c=rep(c("foo", "bar"), each=5))
    r <- DataFrame(a=1:10, b=Rle(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)))
    expect_equal(topdownr:::.colsToRle(d), r)
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
    expect_equal(topdownr:::.droplevels(d[1:5,]), r)
})

test_that(".dropNonInformativeColumns", {
    d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep("foo", 10)))
    r <- DataFrame(a=1:10,
                   c=Rle(rep(c("foo", "bar"), each=5)))
    r2 <- DataFrame(a=1:10,
                    c=Rle(rep(c("foo", "bar"), each=5)),
                    d=Rle(rep("foo", 10)))
    expect_equal(topdownr:::.dropNonInformativeColumns(d), r)
    expect_equal(topdownr:::.dropNonInformativeColumns(d, keep="d"), r2)
})

test_that(".isNumCol", {
    d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep(1:2, each=5)))
    expect_equal(topdownr:::.isNumCol(d), c(TRUE, FALSE, FALSE, TRUE))
})

test_that(".orderByColumns", {
    d <- DataFrame(a=10:1, b=c(1:3, 7:1),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep("foo", 10)))
    expect_error(topdownr:::.orderByColumns(1:10))
    expect_error(topdownr:::.orderByColumns(d, c("foo", "bar")))
    expect_equal(topdownr:::.orderByColumns(d, c("c", "a")), 10:1)
})

test_that(".rbind", {
    l <- rep(c(TRUE, FALSE), 3:2)
    x <- data.frame(a=1:5, b=letters[1:5], c=l, stringsAsFactors=FALSE)
    r <- data.frame(a=1:9, b=letters[1:9], c=c(l, rep(NA, 4)),
                    stringsAsFactors=FALSE)
    y <- list(
        data.frame(a=6:9, b=letters[6:9], stringsAsFactors=FALSE),
        data.frame(a=6:9, c=TRUE, stringsAsFactors=FALSE),
        data.frame(b=letters[6:9], c=TRUE, stringsAsFactors=FALSE),
        data.frame(b=letters[6:9], c=TRUE, d=LETTERS[1:4],
                   stringsAsFactors=FALSE)
    )
    r <- list(
        data.frame(a=1:9, b=letters[1:9], c=c(l, rep(NA, 4)),
                   stringsAsFactors=FALSE),
        data.frame(a=1:9, b=c(letters[1:5], rep(NA_character_, 4)),
                   c=c(l, rep(TRUE, 4)), stringsAsFactors=FALSE),
        data.frame(a=c(1:5, rep(NA_real_, 4)), b=letters[1:9],
                   c=c(l, rep(TRUE, 4)), stringsAsFactors=FALSE),
        data.frame(a=c(1:5, rep(NA_real_, 4)), b=letters[1:9],
                   c=c(l, rep(TRUE, 4)),
                   d=c(rep(NA_character_, 5), LETTERS[1:4]),
                   stringsAsFactors=FALSE)
    )
    expect_error(topdownr:::.rbind(1:10, x))
    expect_error(topdownr:::.rbind(x, 1:10))
    for (i in seq(along=y)) {
        expect_equal(topdownr:::.rbind(x, y[[i]]), r[[i]])
    }
    expect_equal(topdownr:::.rbind(list(x, x, y[[1]])),
                 topdownr:::.rbind(x, x, y[[1]]))
})
