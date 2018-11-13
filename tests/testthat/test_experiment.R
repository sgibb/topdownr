context("experiment")

test_that(".collapseMassList", {
    expect_error(topdownr:::.collapseMassList(1:10))
    expect_error(topdownr:::.collapseMassList(cbind(1:3, 1:3, 1:3)))
    expect_equal(topdownr:::.collapseMassList(cbind(c(10, 20), 1:2)), "10/1 20/2")
})

test_that(".expandMassList", {
    expect_error(topdownr:::.expandMassList(1:10))
    expect_equal(
        topdownr:::.expandMassList("10/1 20/2"),
        matrix(c(10, 20, 1:2), ncol=2, dimnames=list(NULL, c("mz", "z")))
    )
})

test_that(".validMsSettings", {
    expect_error(topdownr:::.validMsSettings(1))
    expect_error(topdownr:::.validMsSettings(FALSE))
    expect_error(topdownr:::.validMsSettings("MS1", family=1))
    expect_error(topdownr:::.validMsSettings("MS1", family="FOO"))
    expect_error(topdownr:::.validMsSettings("MS1", family="Calcium",
                                             version=3))
    expect_error(topdownr:::.validMsSettings("MS1", family="Calcium",
                                             version="foo"))
    expect_true(all(grepl("UVPD.*ActivationTime",
        topdownr:::.validMsSettings("UVPD", family="Calcium",
                                    version="3.2")[, "name"])))
    expect_true(any(grepl("ScanDescription",
        topdownr:::.validMsSettings("MS2", family="Calcium",
                                    version="3.2")[, "name"])))
    expect_false(any(grepl("ScanDescription",
        topdownr:::.validMsSettings("MS2", family="Calcium",
                                    version="3.1")[, "name"])))
})

test_that(".validateMsSetting", {
    expect_true(topdownr:::.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "MS2"
    ))
    expect_true(topdownr:::.validateMsSetting(
        "OrbitrapResolution", c("R15K", "R500K", "R50K"), "ETD"
    ))
    expect_true(topdownr:::.validateMsSetting(
        "FirstMass", 1000, "MS1"
    ))
    expect_true(topdownr:::.validateMsSetting(
        "Microscans", 40L, "MS1"
    ))
    expect_true(topdownr:::.validateMsSetting(
        "MinAgcTarget", TRUE, "MS2"
    ))
    expect_true(grepl("FooBar is not a valid element",
        topdownr:::.validateMsSetting("FooBar", TRUE, "MS2")
    ))
    expect_true(grepl("of type 'MS1'",
        topdownr:::.validateMsSetting("OrbitrapResolution", TRUE, "MS1")
    ))
    expect_true(grepl("could not be 'R40K'",
        topdownr:::.validateMsSetting("OrbitrapResolution", c("R15K", "R40K"), "MS2")
    ))
    expect_true(grepl("should be of class 'integer'",
        topdownr:::.validateMsSetting("Microscans", 10, "MS1")
    ))
})

test_that(".validateMsSettings", {
    expect_true(topdownr:::.validateMsSettings(
        "MS2", list(OrbitrapResolution=c("R15K", "R50K", "R500K"), MinAgcTarget=TRUE)
    ))
    expect_error(topdownr:::.validateMsSettings(
        "MS2", list(FirstMass=1000, OrbitrapResolution=c("R15K", "R50K", "R500K"))
    ), "of type 'MS2'")
})
