context("experiment")

test_that("createExperimentsFragmentOptimisation", {
    ms1 <- data.frame(FirstMass=100, LastMass=200)
    ms2 <- data.frame(OrbitrapResolution="120K", ActivationType="CID",
                      MassList="10/1", stringsAsFactors=FALSE)
    mods <- list(
        MethodModifications=structure(list(
            Modification=structure(list(
                Experiment=structure(
                    list(
                        FullMSScan=list(
                            FirstMass=list(100),
                            LastMass=list(200)
                        ),
                        StartTimeMin=list(0.01),
                        EndTimeMin=list(1)
                    ),
                    ExperimentIndex=0L
                )),
                Order=1L
            ),
            Modification=structure(list(
                    CopyAndAppendExperiment=structure(
                        list(),
                        SourceExperimentIndex=1L
                    ),
                    Experiment=structure(
                        list(),
                        ExperimentIndex=1L
                    )
                ),
                Order=2L
            ),
            Modification=structure(list(
                Experiment=structure(
                    list(
                        TMSnScan=list(
                            MassList=list(
                                MassListRecord=list(
                                    MOverZ=list(10),
                                    Z=list(1)
                                )
                            ),
                            OrbitrapResolution=list("120K"),
                            ActivationType=list("CID"),
                            ScanDescription=list("C1R1")
                        ),
                        StartTimeMin=list(1.01),
                        EndTimeMin=list(2)
                    ),
                    ExperimentIndex=1L
                )),
                Order=3L
            )),
            Version=2L,
            Model="OrbitrapFusionLumos",
            Family="Calcium",
            Type="SL"
        )
    )

    exps <- list("1"=mods, "2"=mods)
    exps[[2]]$MethodModifications[[3]]$Experiment$TMSnScan$ScanDescription  <-
        list("C1R2")
    ## No Start/EndTime
    for (i in seq(along=exps)) {
        for (j in seq(along=exps[[i]]$MethodModifications)) {
            exps[[i]]$MethodModifications[[j]]$Experiment$StartTimeMin <- NULL
            exps[[i]]$MethodModifications[[j]]$Experiment$EndTimeMin <- NULL
        }
    }
    expect_equal(createExperimentsFragmentOptimisation(ms1, ms2,
                    groupBy="replication", randomise=FALSE), exps)
    ## No MassList
    for (i in seq(along=exps)) {
        exps[[i]]$MethodModifications[[3]]$Experiment$TMSnScan$MassList <- NULL
    }
    expect_equal(createExperimentsFragmentOptimisation(ms1,
                    ms2[c("OrbitrapResolution", "ActivationType")],
                    groupBy="replication", randomise=FALSE), exps)

    exps <- list("1"=mods, "2"=mods)
    exps[[2]]$MethodModifications[[3]]$Experiment$TMSnScan$ScanDescription  <-
        list("C1R2")
    expect_equal(
        createExperimentsFragmentOptimisation(
            ms1, ms2, groupBy="replication", scanDuration=1, randomise=FALSE
        ),
        exps
    )
})

test_that(".ms1ConditionToTree", {
    d <- data.frame(FirstMass=100, LastMass=200)
    l1 <- list(Experiment=list(FullMSScan=list(FirstMass=list(100),
                                               LastMass=list(200))))
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    l2$Experiment$StartTimeMin <- list(1)
    l2$Experiment$EndTimeMin <- list(2)
    expect_equal(.ms1ConditionToTree(d, 2, times=NA), l1)
    expect_equal(.ms1ConditionToTree(d, 2, times=1:2), l2)
})

test_that(".copyAndAppendExperiment", {
    l1 <- list(CopyAndAppendExperiment=list(), Experiment=list())
    attr(l1$CopyAndAppendExperiment, "SourceExperimentIndex") <- 0
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    attr(l2$CopyAndAppendExperiment, "SourceExperimentIndex") <- 2
    expect_equal(.copyAndAppendExperiment(2), l1)
    expect_equal(.copyAndAppendExperiment(2, srcId=2), l2)
})

test_that(".tms2ConditionToTree", {
    d <- data.frame(ActivationType="CID", AgcTarget=1e5, stringsAsFactors=FALSE)
    d2 <- cbind(MassList="10/1", d, stringsAsFactors=FALSE)
    d3 <- data.frame(MassList="10/1", stringsAsFactors=FALSE)
    l1 <- list(Experiment=list(TMSnScan=list(ActivationType=list("CID"),
                                             AgcTarget=list(1e5))))
    attr(l1$Experiment, "ExperimentIndex") <- 2
    l2 <- l1
    l2$Experiment$StartTimeMin <- list(1)
    l2$Experiment$EndTimeMin <- list(2)
    l3 <- l1
    l3$Experiment$TMSnScan$MassList <- list(MassListRecord=list(MOverZ=list(10),
                                                       Z=list(1)))
    l4 <- l2
    l4$Experiment$TMSnScan$MassList <- l3$Experiment$TMSnScan$MassList
    l3$Experiment$TMSnScan <- l3$Experiment$TMSnScan[c("MassList",
                                                       "ActivationType",
                                                       "AgcTarget")]
    l4$Experiment$TMSnScan <- l4$Experiment$TMSnScan[c("MassList",
                                                       "ActivationType",
                                                       "AgcTarget")]
    l5 <- list(Experiment=list(TMSnScan=list(MassList=l3$Experiment$TMSnScan$MassList)))
    attr(l5$Experiment, "ExperimentIndex") <- 2

    expect_equal(.tms2ConditionToTree(d, 2, times=NA), l1)
    expect_equal(.tms2ConditionToTree(d, 2, times=1:2), l2)
    expect_equal(.tms2ConditionToTree(d2, 2, times=NA), l3)
    expect_equal(.tms2ConditionToTree(d2, 2, times=1:2), l4)
    expect_equal(.tms2ConditionToTree(d3, 2), l5)
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

test_that("expandTms2Conditions", {
    expect_error(expandTms2Conditions(FOO=1), "FOO is not a valid element")
    expect_error(expandTms2Conditions(ActivationType="CID", AgcTarget=1,
                                     family="FOO"))
    expect_error(expandTms2Conditions(ActivationType="CID", AgcTarget=1,
                                     version="xx"))
    ms2 <- data.frame(MassList="10/1 20/2",
                      ActivationType="CID",
                      AgcTarget=rep(c(1e5, 5e5), 2),
                      OrbitrapResolution=rep(c("R120K", "R60K"), each=2),
                      stringsAsFactors=FALSE)
    expect_equal(expandTms2Conditions(ActivationType="CID",
                                     AgcTarget=c(1e5, 5e5),
                                     OrbitrapResolution=c("R120K", "R60K")),
                 ms2[-1])
    expect_equal(expandTms2Conditions(MassList=cbind(mz=c(10, 20), z=1:2),
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

test_that("validTms2Settings", {
    expect_error(validTms2Settings("FOO"))
    expect_error(validTms2Settings(family=1))
    expect_error(validTms2Settings(family="FOO"))
    expect_error(validTms2Settings(version="xx"))
    m <- validTms2Settings()
    expect_true(class(m) == "matrix")
    expect_equal(colnames(m), c("name", "class", "type"))
    expect_equal(validTms2Settings("All"),
                 validTms2Settings(c("TMS2", "ETD", "CID", "HCD", "UVPD")))
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
        .validMsSettings("TMS2", family="Calcium", version="3.2")[, "name"])))
    expect_false(any(grepl("ScanDescription",
        .validMsSettings("TMS2", family="Calcium", version="3.1")[, "name"])))
})

test_that(".validateMsSetting", {
    expect_true(.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "TMS2"))
    expect_true(.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "ETD"))
    expect_true(.validateMsSetting("FirstMass", 1000, "MS1"))
    expect_true(.validateMsSetting("Microscans", 40L, "MS1"))
    expect_true(.validateMsSetting("MinAgcTarget", TRUE, "TMS2"))
    expect_true(grepl("FooBar is not a valid element",
        .validateMsSetting("FooBar", TRUE, "MS2")))
    expect_true(grepl("of type 'MS1'",
        .validateMsSetting("OrbitrapResolution", TRUE, "MS1")))
    expect_true(grepl("could not be 'R40K'",
        .validateMsSetting("OrbitrapResolution", c("R15K", "R40K"), "TMS2")))
    expect_true(grepl("should be of class 'integer'",
        .validateMsSetting("Microscans", 10, "MS1")))
})

test_that(".validateMsSettings", {
    expect_true(.validateMsSettings("TMS2",
        list(OrbitrapResolution=c("R15K", "R50K", "R500K"), MinAgcTarget=TRUE)))
    expect_error(.validateMsSettings("TMS2",
        list(FirstMass=1000, OrbitrapResolution=c("R15K", "R50K", "R500K"))),
        "of type 'TMS2'"
    )
})
