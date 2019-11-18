context("data.frame")

test_that(".aggregateDataFrame", {
    d <- data.frame(id=1:4, col1=1:4, col2=5:8, col3=LETTERS[1:4],
                    stringsAsFactors=FALSE, row.names=letters[1:4])
    r <- data.frame(id=c(1, 3), col1=c(1.5, 3.5), col2=c(5.5, 7.5),
                    col3=LETTERS[c(1, 3)],
                    stringsAsFactors=FALSE, row.names=letters[c(1, 3)])
    expect_equal(.aggregateDataFrame(d, rep(1:2, each=2),
                                               ignoreNumCols="id"), r)
    expect_equal(.aggregateDataFrame(d, rep(c(10, 2), each=2),
                                               ignoreNumCols="id"), r)
})

test_that(".colsToLogical", {
    d <- DataFrame(a=1:10, b=rep(c("On", "Off"), 5), c=rep(c("foo", "bar"), 5))
    d1 <- DataFrame(
        a=1:10, b=Rle(rep(c("On", "Off"), 5)), c=rep(c("foo", "bar"), 5)
    )
    r <- DataFrame(a=1:10, b=rep(c(TRUE, FALSE), 5), c=rep(c("foo", "bar"), 5))
    expect_equal(.colsToLogical(d), r)
    expect_equal(.colsToLogical(d1), r)
    expect_equal(.colsToLogical(d[, 1, drop=FALSE]), d[, 1, drop=FALSE])
})

test_that(".colsToRle", {
    d <- DataFrame(a=1:10, b=rep(1, 10), c=rep(c("foo", "bar"), each=5))
    r <- DataFrame(a=1:10, b=Rle(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)))
    expect_equal(.colsToRle(d), r)
    expect_equal(.colsToLogical(d[, 1, drop=FALSE]), d[, 1, drop=FALSE])
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
    expect_equal(.droplevels(d[1:5,]), r)
})

test_that(".dropNaColumns", {
    d <- DataFrame(a=1:10, b=NA, c=Rle(rep(c("foo", "bar"), each=5)))
    r <- DataFrame(a=1:10, c=Rle(rep(c("foo", "bar"), each=5)))
    expect_equal(.dropNaColumns(d), r)
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
    expect_equal(.dropNonInformativeColumns(d), r)
    expect_equal(.dropNonInformativeColumns(d, keep="d"), r2)
})

test_that(".isCharacterCol", {
    d <- DataFrame(a=LETTERS[1:10], b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep(1:2, each=5)))
    expect_equal(.isCharacterCol(d), c(TRUE, FALSE, TRUE, FALSE))
})

test_that(".isNumCol", {
    d <- DataFrame(a=1:10, b=factor(rep(1, 10)),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep(1:2, each=5)))
    expect_equal(.isNumCol(d), c(TRUE, FALSE, FALSE, TRUE))
})

test_that(".makeRowNames", {
    d <- data.frame(a=c(1e5, 1e6, 1e7, NA), b=letters[1:4], c=8:11)
    expect_error(.makeRowNames(1:3))
    expect_equal(.makeRowNames(d),
                 c("C1.0e+05_a_08", "C1.0e+06_b_09",
                   "C1.0e+07_c_10", "C0.0e+00_d_11"))
    expect_equal(.makeRowNames(DataFrame(d[, "b", drop=FALSE])),
                 paste0("C", letters[1:4]))
    expect_equal(.makeRowNames(data.frame(a=LETTERS[1:3])),
                 paste0("C", LETTERS[1:3]))
    expect_equal(.makeRowNames(data.frame(a=1:3)),
                 paste0("C", 1:3))
    expect_equal(.makeRowNames(
            data.frame(a=rep(1e5, 4), b=letters[1:4], c=8:11)
        ), c("Ca_08", "Cb_09", "Cc_10", "Cd_11")
    )
    expect_equal(.makeRowNames(data.frame(a=rep(1e5, 10), b="a", c=8)),
                 sprintf("C%02d", 1:10))
})

test_that(".orderByColumns", {
    d <- DataFrame(a=10:1, b=c(1:3, 7:1),
                   c=Rle(rep(c("foo", "bar"), each=5)),
                   d=Rle(rep("foo", 10)))
    expect_error(.orderByColumns(1:10))
    expect_error(.orderByColumns(d, c("foo", "bar")))
    expect_equal(.orderByColumns(d, c("c", "a")), 10:1)
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
                   stringsAsFactors=FALSE),
        DataFrame(a0=6:9, a=6:9, a1=6:9, c=TRUE)
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
                   stringsAsFactors=FALSE),
        DataFrame(a=1:9,
                  b=c(letters[1:5], rep(NA_character_, 4)),
                  c=c(l, rep(TRUE, 4)),
                  a0=c(rep(NA_real_, 5), 6:9),
                  a1=c(rep(NA_real_, 5), 6:9))
    )
    expect_error(.rbind(1:10, x))
    expect_error(.rbind(x, 1:10))
    for (i in seq(along=y)) {
        expect_equal(.rbind(x, y[[i]]), r[[i]])
    }
    expect_equal(.rbind(x), x)
    expect_equal(.rbind(list(x, x, y[[1]])),
                 .rbind(x, x, y[[1]]))
})
