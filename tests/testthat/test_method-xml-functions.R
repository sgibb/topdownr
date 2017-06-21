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
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), groupBy=character(),
                                         replication=1, randomise=FALSE),
               list(data.frame(A=1:2, B="FOO", replication=1,
                               stringsAsFactors=FALSE)), check.attributes=FALSE)
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), groupBy=character(),
                                         randomise=FALSE),
               list(data.frame(A=c(1:2, 1:2), B="FOO", replication=c(1, 1, 2, 2),
                               stringsAsFactors=FALSE)), check.attributes=FALSE)
  set.seed(2017) # set.seed(2017); sample(4) # 4 2 1 3
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), groupBy=character(),
                                         randomise=TRUE),
               list(data.frame(A=c(2, 2, 1, 1), B="FOO", replication=c(2, 1, 1, 2),
                               stringsAsFactors=FALSE)), check.attributes=FALSE)
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), groupBy=c("A"),
                                         randomise=FALSE),
               list("1"=data.frame(A=1, B="FOO", replication=1:2,
                                   stringsAsFactors=FALSE),
                    "2"=data.frame(A=2, B="FOO", replication=1:2,
                                   stringsAsFactors=FALSE)), check.attributes=FALSE)
  expect_equal(topdown:::.ms2Experiments(list(A=1:2, B="FOO"), groupBy=c("A", "replication"),
                                         randomise=FALSE),
               list("1.1"=data.frame(A=1, B="FOO", replication=1,
                                     stringsAsFactors=FALSE),
                    "2.1"=data.frame(A=2, B="FOO", replication=1,
                                     stringsAsFactors=FALSE),
                    "1.2"=data.frame(A=1, B="FOO", replication=2,
                                     stringsAsFactors=FALSE),
                    "2.2"=data.frame(A=2, B="FOO", replication=2,
                                     stringsAsFactors=FALSE)), check.attributes=FALSE)
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

test_that(".startEndTime", {
  r <- data.frame(Type=c("MS1", rep("MS2", 12), "MS1", rep("MS2", 8)),
                  StartTimeMin=seq(0.02, by=0.8, length.out=22),
                  EndTimeMin=seq(0.8, by=0.8, length.out=22),
                  stringsAsFactors=FALSE)
  expect_equal(topdown:::.startEndTime(nMs2=20, nMs2perMs1=12, duration=0.8, gap=0.02),
               r)
  expect_warning(topdown:::.startEndTime(nMs2=201, nMs2perMs1=2),
                 "More than 300 experiments")
})

test_that(".resample", {
  x <- data.frame(A=LETTERS[1:10],
                  B=1:10)
  set.seed(2017) # set.seed(2017); sample(10); # 10  5  4  3  9  8  1  2  6  7
  expect_equal(topdown:::.resample(x), x[c(10, 5:3, 9:8, 1:2, 6:7),])
  expect_equal(topdown:::.resample(x, seq), x)
})

test_that(".xmlHeader", {
  expect_output(topdown:::.xmlHeader(file=""),
                "<\\?xml version=\"1\\.0\" encoding=\"utf-8\"\\?>")
  expect_output(topdown:::.xmlHeader(file="", encoding="utf-16"),
                "<\\?xml version=\"1\\.0\" encoding=\"utf-16\"\\?>")
})

test_that(".xmlTag", {
  expect_output(topdown:::.xmlTag("foo", file=""), "<foo/>")
  expect_output(topdown:::.xmlTag("foo", value="bar", file=""),
                "<foo>bar</foo>")
  expect_output(topdown:::.xmlTag("foo", value="bar", indention=2, file=""),
                "  <foo>bar</foo>")
  expect_output(topdown:::.xmlTag("foo", value="bar", close=FALSE, file=""),
                "<foo>bar")
  expect_output(topdown:::.xmlTag("foo", attrs=c(bar=1), file=""),
                "<foo bar=\"1\"/>")
  expect_output(topdown:::.xmlTag("foo", attrs=c(bar=1), close=FALSE, file=""),
                "<foo bar=\"1\">")
  expect_output(topdown:::.xmlTag("foo", value="bar", attrs=c(x=1, y=2), file=""),
                "<foo x=\"1\" y=\"2\">bar</foo>")
})

test_that(".xmlTagClose", {
  expect_output(topdown:::.xmlTagClose("foo", file=""), "</foo>")
  expect_output(topdown:::.xmlTagClose("foo", indention=2, file=""), "  </foo>")
})

test_that(".xmlListToTags", {
  l <- list(foo="bar", x=1, y=2, z=NA)
  expect_output(topdown:::.xmlListToTags(l, file=""),
                "<foo>bar</foo>\n<x>1</x>\n<y>2</y>")
  expect_output(topdown:::.xmlListToTags(l, indention=2, file=""),
                "  <foo>bar</foo>\n  <x>1</x>\n  <y>2</y>")
  expect_output(topdown:::.xmlListToTags(l, na.rm=FALSE, file=""),
                "<foo>bar</foo>\n<x>1</x>\n<y>2</y>\n<z>NA</z>")
})

test_that(".xmlFullMsScan", {
  expect_output(topdown:::.xmlFullMsScan(list(FirstMass=100, LastMass=200),
                                         file=""),
                paste0("<Modification Order=\"1\">\n",
                       "  <Experiment ExperimentIndex=\"0\">\n",
                       "    <FullMSScan>\n",
                       "      <FirstMass>100</FirstMass>\n",
                       "      <LastMass>200</LastMass>\n",
                       "    </FullMSScan>\n",
                       "  </Experiment>\n",
                       "</Modification>"))
})

test_that(".xmlCopyAndAppendExperiment", {
  expect_output(topdown:::.xmlCopyAndAppendExperiment(order=2, src=1, file=""),
                paste0("<Modification Order=\"2\">\n",
                       "  <CopyAndAppendExperiment SourceExperimentIndex=\"1\"/>\n",
                       "</Modification>"))
})

test_that(".xmlCopyAndAppendExperiment", {
  expect_error(topdown:::.xmlStartEndTime(order=2, times=1, file=""))
  expect_error(topdown:::.xmlStartEndTime(order=2, times=1:3, file=""))
  expect_output(topdown:::.xmlStartEndTime(order=2, times=1:2, file=""),
                paste0("<Modification Order=\"2\">\n",
                       "  <StartTimeMin>1</StartTimeMin>\n",
                       "  <EndTimeMin>2</EndTimeMin>\n",
                       "</Modification>"))
})
