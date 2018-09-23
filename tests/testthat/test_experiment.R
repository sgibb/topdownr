context("experiment")

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
