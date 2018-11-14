context("experiment")

test_that(".ms1ConditionToTree", {
    d <- data.frame(FirstMass=100, LastMass=200)
    l1 <- list(Experiment=list(FullMSScan=list(FirstMass=list(100),
                                               LastMass=list(200))))
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    l2$Experiment$StartTimeMin <- list(1)
    l2$Experiment$EndTimeMin <- list(2)
    expect_equal(.ms1ConditionToTree(d, 2, times=NULL), l1)
    expect_equal(.ms1ConditionToTree(d, 2, times=1:2), l2)
})

test_that(".ms1CopyAndAppendExperiment", {
    l1 <- list(CopyAndAppendExperiment=list(), Experiment=list())
    attr(l1$CopyAndAppendExperiment, "SourceExperimentIndex") <- 0
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    l2$Experiment$StartTimeMin <- list(1)
    l2$Experiment$EndTimeMin <- list(2)
    expect_equal(.ms1CopyAndAppendExperiment(2, times=NULL), l1)
    expect_equal(.ms1CopyAndAppendExperiment(2, times=1:2), l2)
})

test_that(".ms2ConditionToTree", {
    d <- data.frame(ActivationType="CID", AgcTarget=1e5, stringsAsFactors=FALSE)
    d2 <- cbind(MassList="10/1", d, stringsAsFactors=FALSE)
    l1 <- list(Experiment=list(TMSnScan=list(ActivationType=list("CID"),
                                             AgcTarget=list(1e5))))
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    l2$Experiment$StartTimeMin <- list(1)
    l2$Experiment$EndTimeMin <- list(2)
    l3 <- l1
    l3$Experiment$MassList <- list(MassListRecord=list(MOverZ=list(10),
                                                       Z=list(1)))
    l4 <- l2
    l4$Experiment$MassList <- l3$Experiment$MassList
    expect_equal(.ms2ConditionToTree(d, 2, times=NULL), l1)
    expect_equal(.ms2ConditionToTree(d, 2, times=1:2), l2)
    expect_equal(.ms2ConditionToTree(d2, 2, times=NULL), l3)
    expect_equal(.ms2ConditionToTree(d2, 2, times=1:2), l4)
})

test_that(".collapseMassList", {
    ml1 <- list(MassListRecord=list(MOverZ=list(10), Z=list(1)))
    ml2 <- list(
        MassListRecord=list(MOverZ=list(10), Z=list(1)),
        MassListRecord=list(MOverZ=list(20), Z=list(2))
    )
    expect_equal(.massListToTree("10/1"), ml1)
    expect_equal(.massListToTree("10/1 20/2"), ml2)
})

test_that(".collapseMassList", {
    expect_error(.collapseMassList(1:10))
    expect_error(.collapseMassList(cbind(1:3, 1:3, 1:3)))
    expect_equal(.collapseMassList(cbind(c(10, 20), 1:2)), "10/1 20/2")
})

test_that(".expandMassList", {
    expect_error(.expandMassList(1:10))
    expect_equal(
        .expandMassList("10/1 20/2"),
        matrix(c(10, 20, 1:2), ncol=2, dimnames=list(NULL, c("mz", "z")))
    )
})

test_that("expandMs1Conditions", {
    expect_error(expandMs1Conditions(FOO=1), "FOO is not a valid element")
    expect_error(expandMs1Conditions(FirstMass=1, family="FOO"))
    expect_error(expandMs1Conditions(FirstMass=1, version="xx"))
    expect_equal(expandMs1Conditions(FirstMass=100, LastMass=200),
                 data.frame(FirstMass=100, LastMass=200), )
})

test_that("expandMs2Conditions", {
    expect_error(expandMs2Conditions(FOO=1), "FOO is not a valid element")
    expect_error(expandMs2Conditions(ActivationType="CID", AgcTarget=1,
                                     family="FOO"))
    expect_error(expandMs2Conditions(ActivationType="CID", AgcTarget=1,
                                     version="xx"))
    ms2 <- data.frame(MassList="10/1 20/2",
                      ActivationType="CID",
                      AgcTarget=rep(c(1e5, 5e5), 2),
                      OrbitrapResolution=rep(c("R120K", "R60K"), each=2),
                      stringsAsFactors=FALSE)
    expect_equal(expandMs2Conditions(ActivationType="CID",
                                     AgcTarget=c(1e5, 5e5),
                                     OrbitrapResolution=c("R120K", "R60K")),
                 ms2[-1])
    expect_equal(expandMs2Conditions(MassList=cbind(mz=c(10, 20), z=1:2),
                                     ActivationType="CID",
                                     AgcTarget=c(1e5, 5e5),
                                     OrbitrapResolution=c("R120K", "R60K")),
                 ms2)
})

test_that("validMs1Settings", {
    expect_error(validMs1Settings(family=1))
    expect_error(validMs1Settings(family="FOO"))
    expect_error(validMs1Settings(version="xx"))
    m <- validMs1Settings()
    expect_true(class(m) == "matrix")
    expect_equal(colnames(m), c("name", "class", "type"))
    expect_equal(m[1:3, "name"], c("FirstMass", "LastMass", "Microscans"))
    expect_equal(m[1:3, "class"], c("double", "double", "integer"))
    expect_equal(m[1:3, "type"], rep("MS1", 3))
})

test_that("validMs2Settings", {
    expect_error(validMs2Settings("FOO"))
    expect_error(validMs2Settings(family=1))
    expect_error(validMs2Settings(family="FOO"))
    expect_error(validMs2Settings(version="xx"))
    m <- validMs2Settings()
    expect_true(class(m) == "matrix")
    expect_equal(colnames(m), c("name", "class", "type"))
    expect_equal(validMs2Settings("All"),
                 validMs2Settings(c("MS2", "ETD", "CID", "HCD", "UVPD")))
})

test_that(".validMsSettings", {
    expect_error(.validMsSettings(1))
    expect_error(.validMsSettings(FALSE))
    expect_error(.validMsSettings("MS1", family=1))
    expect_error(.validMsSettings("MS1", family="FOO"))
    expect_error(.validMsSettings("MS1", family="Calcium",
                                             version=3))
    expect_error(.validMsSettings("MS1", family="Calcium", version="foo"))
    expect_true(all(grepl("UVPD.*ActivationTime",
        .validMsSettings("UVPD", family="Calcium", version="3.2")[, "name"])))
    expect_true(any(grepl("ScanDescription",
        .validMsSettings("MS2", family="Calcium", version="3.2")[, "name"])))
    expect_false(any(grepl("ScanDescription",
        .validMsSettings("MS2", family="Calcium", version="3.1")[, "name"])))
})

test_that(".validateMsSetting", {
    expect_true(.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "MS2"))
    expect_true(.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "ETD"))
    expect_true(.validateMsSetting("FirstMass", 1000, "MS1"))
    expect_true(.validateMsSetting("Microscans", 40L, "MS1"))
    expect_true(.validateMsSetting("MinAgcTarget", TRUE, "MS2"))
    expect_true(grepl("FooBar is not a valid element",
        .validateMsSetting("FooBar", TRUE, "MS2")))
    expect_true(grepl("of type 'MS1'",
        .validateMsSetting("OrbitrapResolution", TRUE, "MS1")))
    expect_true(grepl("could not be 'R40K'",
        .validateMsSetting("OrbitrapResolution", c("R15K", "R40K"), "MS2")))
    expect_true(grepl("should be of class 'integer'",
        .validateMsSetting("Microscans", 10, "MS1")))
})

test_that(".validateMsSettings", {
    expect_true(.validateMsSettings("MS2",
        list(OrbitrapResolution=c("R15K", "R50K", "R500K"), MinAgcTarget=TRUE)))
    expect_error(.validateMsSettings("MS2",
        list(FirstMass=1000, OrbitrapResolution=c("R15K", "R50K", "R500K"))),
        "of type 'MS2'"
    )
})
