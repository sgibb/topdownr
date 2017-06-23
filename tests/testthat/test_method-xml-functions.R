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
               list("1:1"=data.frame(A=1, B="FOO", replication=1,
                                     stringsAsFactors=FALSE),
                    "1:2"=data.frame(A=1, B="FOO", replication=2,
                                     stringsAsFactors=FALSE),
                    "2:1"=data.frame(A=2, B="FOO", replication=1,
                                     stringsAsFactors=FALSE),
                    "2:2"=data.frame(A=2, B="FOO", replication=2,
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
  x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), na=rep(c(1, NA), 4),
                  stringsAsFactors=FALSE)
  expect_equal(topdown:::.groupExperimentsBy(x, "LE"), split(x, x$LE))
  expect_equal(topdown:::.groupExperimentsBy(x, c("ID", "LE")),
               split(x, interaction(as.list(x[, c("ID", "LE")]),
                                    sep=":", lex.order=TRUE)))
  expect_equal(topdown:::.groupExperimentsBy(x, c("ID", "na")),
               setNames(split(x, x$ID), c("1:1", "2:NA")))
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

test_that(".xmlTMSnScan", {
  expect_output(topdown:::.xmlTMSnScan(data.frame(OrbitrapResolution="R120K",
                                                  ActivationType="CID",
                                                  CollisionEnergy=5,
                                                  stringsAsFactors=FALSE),
                                       mz=c(100, 120),
                                       z=2:3,
                                       order=3,
                                       idx=2,
                                       file=""),
                paste0("<Modification Order=\"3\">\n",
                       "  <Experiment ExperimentIndex=\"2\">\n",
                       "    <TMSnScan>\n",
                       "      <OrbitrapResolution>R120K</OrbitrapResolution>\n",
                       "      <ActivationType>CID</ActivationType>\n",
                       "      <MassList CollisionEnergyCID=\"true\">\n",
                       "        <MassListRecord>\n",
                       "          <MOverZ>100</MOverZ>\n",
                       "          <Z>2</Z>\n",
                       "          <CollisionEnergyCID>5</CollisionEnergyCID>\n",
                       "        </MassListRecord>\n",
                       "        <MassListRecord>\n",
                       "          <MOverZ>120</MOverZ>\n",
                       "          <Z>3</Z>\n",
                       "          <CollisionEnergyCID>5</CollisionEnergyCID>\n",
                       "        </MassListRecord>\n",
                       "      </MassList>\n",
                       "    </TMSnScan>\n",
                       "  </Experiment>\n",
                       "</Modification>"))
})

test_that(".xmlCopyAndAppendExperiment", {
  expect_output(topdown:::.xmlCopyAndAppendExperiment(order=2, src=1, file=""),
                paste0("<Modification Order=\"2\">\n",
                       "  <CopyAndAppendExperiment SourceExperimentIndex=\"1\"/>\n",
                       "</Modification>"))
})

test_that(".xmlStartEndTime", {
  expect_error(topdown:::.xmlStartEndTime(order=2, times=1, file=""))
  expect_error(topdown:::.xmlStartEndTime(order=2, times=1:3, file=""))
  expect_output(topdown:::.xmlStartEndTime(order=2, times=1:2, file=""),
                paste0("<Modification Order=\"2\">\n",
                       "  <StartTimeMin>1</StartTimeMin>\n",
                       "  <EndTimeMin>2</EndTimeMin>\n",
                       "</Modification>"))
})

test_that(".xmlMassList", {
  expect_output(topdown:::.xmlMassList(c(100, 120), 2:3, file=""),
                paste0("<MassList>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>100</MOverZ>\n",
                       "    <Z>2</Z>\n",
                       "  </MassListRecord>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>120</MOverZ>\n",
                       "    <Z>3</Z>\n",
                       "  </MassListRecord>\n",
                       "</MassList>"))
  expect_output(topdown:::.xmlMassList(c(100, 120), 14:15, file=""),
                paste0("<MassList>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>100</MOverZ>\n",
                       "    <Z>10</Z>\n",
                       "  </MassListRecord>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>120</MOverZ>\n",
                       "    <Z>10</Z>\n",
                       "  </MassListRecord>\n",
                       "</MassList>"))
  expect_output(topdown:::.xmlMassList(c(100, 120, 140), 2:4,
                                       indention=2, file=""),
                paste0("  <MassList>\n",
                       "    <MassListRecord>\n",
                       "      <MOverZ>100</MOverZ>\n",
                       "      <Z>2</Z>\n",
                       "    </MassListRecord>\n",
                       "    <MassListRecord>\n",
                       "      <MOverZ>120</MOverZ>\n",
                       "      <Z>3</Z>\n",
                       "    </MassListRecord>\n",
                       "    <MassListRecord>\n",
                       "      <MOverZ>140</MOverZ>\n",
                       "      <Z>4</Z>\n",
                       "    </MassListRecord>\n",
                       "  </MassList>"))
  expect_output(topdown:::.xmlMassList(c(100, 120), 2:3, 30, type="CID", file=""),
                paste0("<MassList CollisionEnergyCID=\"true\">\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>100</MOverZ>\n",
                       "    <Z>2</Z>\n",
                       "    <CollisionEnergyCID>30</CollisionEnergyCID>\n",
                       "  </MassListRecord>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>120</MOverZ>\n",
                       "    <Z>3</Z>\n",
                       "    <CollisionEnergyCID>30</CollisionEnergyCID>\n",
                       "  </MassListRecord>\n",
                       "</MassList>"))
  expect_output(topdown:::.xmlMassList(c(100, 120), 2:3, 50, type="HCD", file=""),
                paste0("<MassList CollisionEnergyHCD=\"true\">\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>100</MOverZ>\n",
                       "    <Z>2</Z>\n",
                       "    <CollisionEnergyHCD>50</CollisionEnergyHCD>\n",
                       "  </MassListRecord>\n",
                       "  <MassListRecord>\n",
                       "    <MOverZ>120</MOverZ>\n",
                       "    <Z>3</Z>\n",
                       "    <CollisionEnergyHCD>50</CollisionEnergyHCD>\n",
                       "  </MassListRecord>\n",
                       "</MassList>"))
})

test_that(".xmlMassListRecord", {
  expect_output(topdown:::.xmlMassListRecord(100, 3, file=""),
                paste0("<MassListRecord>\n",
                       "  <MOverZ>100</MOverZ>\n",
                       "  <Z>3</Z>\n",
                       "</MassListRecord>"))
  expect_output(topdown:::.xmlMassListRecord(100, 13, file=""),
                paste0("<MassListRecord>\n",
                       "  <MOverZ>100</MOverZ>\n",
                       "  <Z>10</Z>\n",
                       "</MassListRecord>"))
  expect_output(topdown:::.xmlMassListRecord(100, 3, indention=2, file=""),
                paste0("  <MassListRecord>\n",
                       "    <MOverZ>100</MOverZ>\n",
                       "    <Z>3</Z>\n",
                       "  </MassListRecord>"))
  expect_output(topdown:::.xmlMassListRecord(20, 5, energy=30, type="CID",
                                             file=""),
                paste0("<MassListRecord>\n",
                       "  <MOverZ>20</MOverZ>\n",
                       "  <Z>5</Z>\n",
                       "  <CollisionEnergyCID>30</CollisionEnergyCID>\n",
                       "</MassListRecord>"))
  expect_output(topdown:::.xmlMassListRecord(20, 5, energy=50, type="HCD",
                                             file=""),
                paste0("<MassListRecord>\n",
                       "  <MOverZ>20</MOverZ>\n",
                       "  <Z>5</Z>\n",
                       "  <CollisionEnergyHCD>50</CollisionEnergyHCD>\n",
                       "</MassListRecord>"))
})

test_that("writeMethodXmls", {
  expect_error(topdown:::writeMethodXmls(list(1:3)),
               ".*ms1Settings.* has to be a named list")
  expect_error(topdown:::writeMethodXmls(list(foo=1:3)),
               "is/are no valid MS1 tag")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(), list(1:3)),
               ".*ms2Settings.* has to be a named list")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(), list(foo=1:3)),
               "is/are no valid MS2 tag")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(),
                                         defaultMs2Settings(),
                                         groupBy="foo"),
               "Items of .*groupBy.* have to be one or more of:")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(),
                                         defaultMs2Settings(),
                                         mz=1),
               ".*mz.* has to be a matrix")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(),
                                         defaultMs2Settings(),
                                         mz=cbind(1:3, 1:3, 1:3)),
               ".*mz.* has to be a matrix with two columns")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(),
                                         defaultMs2Settings(),
                                         mz=matrix(nrow=0, ncol=2)),
               ".*mz.* has to have at least one row")
  expect_error(topdown:::writeMethodXmls(defaultMs1Settings(),
                                         defaultMs2Settings(),
                                         mz=cbind(1:3, 1:3),
                                         pattern="foo.xml"),
               " has to contain '%s' to be replaced")

  xml <- c('<?xml version="1.0" encoding="utf-8"?>',
           '<MethodModifications Version="1" Model="OrbitrapFusion" Family="Calcium" Type="SL">',
           '<Modification Order="1">',
           '  <Experiment ExperimentIndex="0">',
           '    <FullMSScan>',
           '      <FirstMass>100</FirstMass>',
           '    </FullMSScan>',
           '  </Experiment>',
           '</Modification>',
           '<Modification Order="2">',
           '  <CopyAndAppendExperiment SourceExperimentIndex="1"/>',
           '</Modification>',
           '<Modification Order="3">',
           '  <StartTimeMin>0.01</StartTimeMin>',
           '  <EndTimeMin>0.5</EndTimeMin>',
           '</Modification>',
           '<Modification Order="4">',
           '  <StartTimeMin>0.51</StartTimeMin>',
           '  <EndTimeMin>1</EndTimeMin>',
           '</Modification>',
           '<Modification Order="5">',
           '  <StartTimeMin>1.01</StartTimeMin>',
           '  <EndTimeMin>1.5</EndTimeMin>',
           '</Modification>',
           '<Modification Order="6">',
           '  <Experiment ExperimentIndex="1">',
           '    <TMSnScan>',
           '      <ActivationType>ETD</ActivationType>',
           '      <AgcTarget>10000</AgcTarget>',
           '      <ETDReactionTime>10</ETDReactionTime>',
           '      <replication>1</replication>',
           '      <MassList>',
           '        <MassListRecord>',
           '          <MOverZ>100.0001</MOverZ>',
           '          <Z>2</Z>',
           '        </MassListRecord>',
           '      </MassList>',
           '    </TMSnScan>',
           '  </Experiment>',
           '</Modification>',
           '<Modification Order="7">',
           '  <Experiment ExperimentIndex="2">',
           '    <TMSnScan>',
           '      <ActivationType>ETD</ActivationType>',
           '      <AgcTarget>20000</AgcTarget>',
           '      <ETDReactionTime>10</ETDReactionTime>',
           '      <replication>1</replication>',
           '      <MassList>',
           '        <MassListRecord>',
           '          <MOverZ>100.0002</MOverZ>',
           '          <Z>2</Z>',
           '        </MassListRecord>',
           '      </MassList>',
           '    </TMSnScan>',
           '  </Experiment>',
           '</Modification>',
           '</MethodModifications>')

  tdir <- tempdir()
  writeMethodXmls(list(FirstMass=100),
                  list(ActivationType="ETD",
                       AgcTarget=c(10000, 20000),
                       ETDReactionTime=c(10, 20)),
                  groupBy="ETDReactionTime",
                  replications=1,
                  mz=cbind(100, 2),
                  randomise=FALSE,
                  pattern=file.path(tdir, "method_%s.xml"))

  expect_equal(readLines(file.path(tdir, "method_10.xml")), xml)
  expect_equal(readLines(file.path(tdir, "method_20.xml")),
               gsub("<ETDReactionTime>10", "<ETDReactionTime>20", xml))
  unlink(list.files(tdir, pattern="method.*\\.xml$", full.names=TRUE))
})
