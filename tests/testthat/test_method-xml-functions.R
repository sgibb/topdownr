context("method-xml")

test_that(".massLabel", {
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245)),
               c(750.0001, 1000.8245))
  expect_equal(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e5),
               c(750.00001, 1000.80245))
  expect_error(topdown:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e3),
               "at least two digits more than")
})

test_that(".ms2Experiments", {
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO")),
               data.frame(A=c(1:2, 1:2), B="FOO", replication=c(1, 1, 2, 2),
                          stringsAsFactors=FALSE), check.attributes=FALSE)
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), replication=1),
               data.frame(A=1:2, B="FOO", replication=1,
                          stringsAsFactors=FALSE), check.attributes=FALSE)
})

test_that(".replaceZeroETDReactionTime", {
  x <- data.frame(ActivationType="ETD",
                  ETDReactionTime=c(0, 0, 0, 1:3),
                  ETDSupplementalActivation="ETcid",
                  ETDSupplementalActivationEnergy=1:6,
                  stringsAsFactors=FALSE)
  r <- data.frame(ActivationType=rep(c("CID","ETD"), each=3),
                  ETDReactionTime=c(rep(NA, 3), 1:3),
                  ETDSupplementalActivation=rep(c(NA, "ETcid"), each=3),
                  ETDSupplementalActivationEnergy=c(rep(NA, 3), 4:6),
                  CollisionEnergy=c(1:3, rep(NA, 3)),
                  stringsAsFactors=FALSE)

  expect_equal(topdown:::.replaceZeroETDReactionTime(x), r)
})

test_that(".groupExperimentsBy", {
  x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), stringsAsFactors=FALSE)
  expect_equal(topdown:::.groupExperimentsBy(x, "LE"), split(x, x$LE))
  expect_equal(topdown:::.groupExperimentsBy(x, c("ID", "LE")),
               split(x, interaction(as.list(x[, c("ID", "LE")]))))
})
