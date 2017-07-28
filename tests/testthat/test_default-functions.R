context("settings-functions")

test_that(".defaultSettings", {
  expect_equal(topdown:::.defaultSettings(default=list(foo=2, bar="x"),
                                          valid=c("foo", "bar")),
               list(foo=2, bar="x"))
  expect_equal(topdown:::.defaultSettings(foo=1,
                                          default=list(foo=2, bar="x"),
                                          valid=c("foo", "bar")),
               list(foo=1, bar="x"))
  expect_error(topdown:::.defaultSettings(foo=1, BAR=2, default=list(foo=2, bar="x"),
                                          valid=c("foo", "bar")),
               "The following setting\\(s\\) is/are not valid: BAR")
})

test_that("defaultMs1Settings", {
  expect_equal(names(defaultMs1Settings()),
               c("FirstMass", "LastMass", "Microscans"))
})

test_that("defaultMs2Settings", {
  expect_equal(names(defaultMs2Settings()),
               c("ActivationType", "AgcTarget", "ETDReagentTarget",
                 "ETDReactionTime", "ETDSupplementalActivation",
                 "ETDSupplementalActivationEnergy", "IsolationWindow",
                 "MaxITTimeInMS", "Microscans", "OrbitrapResolution"))
})
