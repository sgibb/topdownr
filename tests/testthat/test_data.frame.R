context("data.frame")

test_that(".aggregateDataFrame", {
  d <- data.frame(id=1:4, col1=1:4, col2=5:8, col3=LETTERS[1:4],
                  stringsAsFactors=FALSE, row.names=letters[1:4])
  r <- data.frame(col1=c(1.5, 3.5), col2=c(5.5, 7.5), col3=LETTERS[c(1, 3)],
                  stringsAsFactors=FALSE, row.names=letters[c(1, 3)])
  expect_equal(topdown:::.aggregateDataFrame(d, rep(1:2, each=2),
                                             ignoreCols="id"), r)
})
