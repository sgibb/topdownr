context("utils")

test_that(".allIdentical", {
    expect_true(topdownr:::.allIdentical(rep(1L, 3)))
    expect_false(topdownr:::.allIdentical(1:3))
})

test_that(".camelCase", {
    expect_equal(topdownr:::.camelCase(c("Monoisotopic M/Z", "SPS Mass 2",
                                        "RT (min)", "MSLevel", "peaksCount",
                                        "Multi.Inject.Info", "RF.Comp...ppm",
                                        "SupplementalActivationCE",
                                        "TIC", "UseCalibratedUVPDTime",
                                        "UseCalibratedUVPDTimeMS2", "EThcD")),
                 c("MonoisotopicMz", "SpsMass2", "RtMin", "MsLevel",
                   "PeaksCount", "MultiInjectInfo", "RfCompPpm",
                   "SupplementalActivationCe", "Tic",
                   "UseCalibratedUvpdTime", "UseCalibratedUvpdTimeMs2",
                   "Ethcd"))
})

test_that("cat0", {
    expect_output(topdownr:::cat0("foo", "bar"), "foobar")
})

test_that("characterToLogical", {
    expect_error(topdownr:::.characterToLogical(1:2))
    expect_equal(topdownr:::.characterToLogical(c("FOO", "BAR")),
                 c("FOO", "BAR"))
    expect_equal(topdownr:::.characterToLogical(c("On", "Off", "NA")),
                 c(TRUE, FALSE, NA))
    expect_equal(topdownr:::.characterToLogical(c("On", "Off", "NA"),
                                                na.strings="FOO"),
                 c("On", "Off", "NA"))
    expect_equal(topdownr:::.characterToLogical(c("true", "False", "N/A")),
                 c(TRUE, FALSE, NA))
    expect_equal(topdownr:::.characterToLogical(c("true", "False", "N/A",
                                                  NA_character_)),
                 c(TRUE, FALSE, NA, NA))
    expect_equal(topdownr:::.characterToLogical(c("true", "on", "ON")),
                 rep(TRUE, 3))
    expect_equal(topdownr:::.characterToLogical(Rle(c("true", "False", "N/A"))),
                 c(TRUE, FALSE, NA))
    expect_equal(topdownr:::.characterToLogical(Rle(c("true", "on", "ON"))),
                 rep(TRUE, 3))
})

test_that(".fileExt", {
    f <- c("foo.bar", "foo.bar.gz")
    expect_equal(topdownr:::.fileExt(f), c("bar", "bar"))
    expect_equal(topdownr:::.fileExt(f, compression=FALSE), c("bar", "gz"))
})

test_that(".filterStringToId", {
    expect_error(topdownr:::.filterStringToId(1:3))
    expect_equal(topdownr:::.filterStringToId(
    c("FTMS + p NSI sa Full ms2 560.6219@etd50.00@cid7.00 [160.0000-2000.0000]",
      "FTMS + p NSI Full ms2 560.6010@hcd35.00 [160.0000-2000.0000]")),
      c(219L, 10L))
})

test_that(".fixFilterString", {
    fs <- c("FTMS + p NSI Full ms2 162.0004@cid28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 162.0004@hcd28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 162.0004@hcd28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 162.0006@cid35.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 1162.0008@cid28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 1162.0010@hcd28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 1162.0010@hcd28.00 [100.0000-2000.0000]",
            "FTMS + p NSI Full ms2 1162.0010@cid35.00 [100.0000-2000.0000]")
    ffs <- fs
    ffs[2:3] <- "FTMS + p NSI Full ms2 162.0005@hcd28.00 [100.0000-2000.0000]"
    ffs[6:7] <- "FTMS + p NSI Full ms2 1162.0009@hcd28.00 [100.0000-2000.0000]"
    expect_equal(topdownr:::.fixFilterString(fs), ffs)
})

test_that(".fixFilterStringId", {
    expect_equal(topdownr:::.fixFilterStringId(c(1, 2, 2, 4)), 1:4)
    expect_equal(topdownr:::.fixFilterStringId(c(1, 3, 3, 4)), 1:4)
    expect_equal(topdownr:::.fixFilterStringId(c(4, 4, 6, 8, 10, 10)),
                 c(4:6, 8:10))
    expect_equal(topdownr:::.fixFilterStringId(c(5, 5, 6, 8, 9, 9)),
                 c(4:6, 8:10))
})

test_that(".formatNumbers", {
    expect_equal(topdownr:::.formatNumbers(1:10), sprintf("%02d", 1:10))
    expect_equal(topdownr:::.formatNumbers(c(1, 100.1)),
                 sprintf("%06.2f", c(1, 100.1)))
    expect_equal(topdownr:::.formatNumbers(c(1, 100) + 0.1),
                 sprintf("%06.2f", c(1, 100) + 0.1))
    expect_equal(topdownr:::.formatNumbers(c(1, 100) + 0.1, asInteger=TRUE),
                 sprintf("%03d", c(1L, 100L)))
    expect_equal(topdownr:::.formatNumbers(c(1, 100) + 0.1, asInteger=FALSE),
                 sprintf("%06.2f", c(1, 100) + 0.1))
    expect_equal(topdownr:::.formatNumbers(c(1, 100) + 0.1, asInteger=NA),
                 sprintf("%06.2f", c(1, 100) + 0.1))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6) + 0.1),
                 sprintf("%010.2f", c(1, 1000, 1e6) + 0.1))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6)),
                 sprintf("%.1e", c(1, 1000, 1e6)))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6), asInteger=TRUE),
                 sprintf("%.1e", c(1, 1000, 1e6)))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6), nScientific=10),
                 sprintf("%07d", c(1, 1000, 1e6)))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6, NA)),
                 sprintf("%.1e", c(1, 1000, 1e6, NA)))
    expect_equal(topdownr:::.formatNumbers(c(1, 1000, 1e6, NA), na2zero=TRUE),
                 sprintf("%.1e", c(1, 1000, 1e6, 0)))
    expect_equal(topdownr:::.formatNumbers(c(1, -1000)), c("00001", "-1000"))
    expect_equal(topdownr:::.formatNumbers(Rle(c(1.1, 2.1), 2:3)),
                 rep(c("1.10", "2.10"), 2:3))
})

test_that(".fragmentationMethod", {
    d <- expand.grid(EtdActivation=c(NA, 1),
                     CidActivation=c(NA, 1),
                     HcdActivation=c(NA, 1),
                     UvpdActivation=c(NA, 1))
    expect_error(topdownr:::.fragmentationMethod(cbind(d, foo=1L)))
    expect_equal(topdownr:::.fragmentationMethod(d),
                 c("None", "ETD", "CID", "ETcid", "HCD", "EThcd", "CID/HCD",
                   "ETD/CID/HCD", "UVPD", rep(NA_character_, 6), "All"))
})

test_that(".groupBy", {
    x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), na=rep(c(1, NA), 4),
                    stringsAsFactors=FALSE)
    expect_error(topdownr:::.groupBy(1:10, "LE"))
    expect_equal(topdownr:::.groupBy(x, "LE"), split(x, x$LE))
    expect_equal(topdownr:::.groupBy(x, c("ID", "LE")),
                 split(x, interaction(as.list(x[, c("ID", "LE")]),
                                      sep=":", lex.order=TRUE)))
    expect_equal(topdownr:::.groupBy(x, c("ID", "na")),
                 setNames(split(x, x$ID), c("1:1", "2:NA")))
})

test_that(".groupByLabels", {
    x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=2), na=rep(c(1, NA), 4),
                    stringsAsFactors=FALSE)
    DF <- DataFrame(ID=1:2, LE=rep(LETTERS[1:4], each=2))
    expect_error(topdownr:::.groupByLabels(1:10, "LE"), "valid column names")
    expect_equal(topdownr:::.groupByLabels(x, "LE"), x$LE)
    expect_equal(topdownr:::.groupByLabels(x, c("ID", "LE")),
                 paste(1:2, rep(LETTERS[1:4], each=2), sep=":"))
    expect_equal(topdownr:::.groupByLabels(DF),
                 paste(1:2, rep(LETTERS[1:4], each=2), sep=":"))
    expect_equal(topdownr:::.groupByLabels(DF, c("ID", "LE")),
                 paste(1:2, rep(LETTERS[1:4], each=2), sep=":"))
    expect_equal(topdownr:::.groupByLabels(x, c("ID", "na")),
                 paste(rep(1:2, 4), rep(c(1, NA), 4), sep=":"))
    expect_equal(topdownr:::.groupByLabels(x, c("ID", "na"), sep="_"),
                 paste(rep(1:2, 4), rep(c(1, NA), 4), sep="_"))
})

test_that(".groupId", {
    x <- data.frame(ID=1:2, LE=rep(LETTERS[1:4], each=4),
                    stringsAsFactors=FALSE)
    expect_error(topdownr:::.groupId(1:10, "LE"))
    expect_equal(topdownr:::.groupId(x, "LE"), rep(1:4, each=4))
    expect_equal(topdownr:::.groupId(x, c("ID", "LE")),
                 rep(1:2, 8) + rep(seq(0, 6, by=2), each=4))
})

test_that(".hft", {
    expect_equal(topdownr:::.hft(letters[1:6]), letters[1:6])
    expect_equal(topdownr:::.hft(letters[1:26]),
                 c("a", "b", "c", "...", "x", "y", "z"))
    expect_equal(topdownr:::.hft(letters[1:26], fill=NULL, n=4),
                 c("a", "b", "c", "d", "w", "x", "y", "z"))
})

test_that(".logmsg", {
    rxDate <-paste("^\\[20[0-9]{2}-[01][0-9]-[0-3][0-9]",
                   "[0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\]")

    expect_true(grepl(paste(rxDate, "foo$"), topdownr:::.logmsg("foo")))
    expect_true(grepl(paste(rxDate, "foobar$"),
                      topdownr:::.logmsg("foo", "bar")))
})

test_that(".makeNames", {
    x <- rep(LETTERS[1:3], c(2, 1, 10))
    expect_equal(topdownr:::.makeNames(x),
                 c("A:1", "A:2", "B", sprintf("C:%02d", 1:10)))
    expect_equal(topdownr:::.makeNames(x, sep="_", prefix="D"),
                 c("DA_1", "DA_2", "DB", sprintf("DC_%02d", 1:10)))
})

test_that(".massLabel", {
    expect_equal(topdownr:::.massLabel(c(750, 1000.76), c(1, 245)),
                 c(750.0001, 1000.8245))
    expect_equal(topdownr:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e5),
                 c(750.00001, 1000.80245))
    expect_equal(topdownr:::.massLabel(1, 1:999),
                 as.double(sprintf("1.%04d", 1:999)))
    expect_error(topdownr:::.massLabel(c(750, 1000.76), c(1, 245), divisor=1e3),
                 "at least two digits more than")
})

test_that(".massLabelToId", {
    expect_equal(topdownr:::.massLabelToId(c("750.0001", "1000.8245")),
                 c(1, 245))
    expect_equal(topdownr:::.massLabelToId(c("750.001", "1000.824"), 2),
                 c(1, 24))
    expect_equal(topdownr:::.massLabelToId(sprintf("1000.%04d", 1:999)),
                 c(1:999))
})

test_that(".msg", {
    expect_message(topdownr:::.msg(TRUE, "foobar"), "foobar")
    expect_message(topdownr:::.msg(TRUE, "foo", "bar"), "foobar")
    expect_silent(topdownr:::.msg(FALSE, "foobar"))
})

test_that(".ndigits", {
    expect_equal(topdownr:::.ndigits(rep(10^(1:6), each=2) - c(0, 1)),
                 rep(2:7, each=2) - c(0, 1))
    expect_equal(topdownr:::.ndigits(-c(1, 10)), 1:2)
    expect_equal(topdownr:::.ndigits(0), 1)
    expect_equal(topdownr:::.ndigits(c(NA, 30)), 2)
    expect_equal(topdownr:::.ndigits(c(NA)), 1)
})

test_that(".nrows", {
    expect_error(topdownr:::.nrows(matrix(nrow=2, ncol=2)))
    expect_equal(topdownr:::.nrows(list(matrix(nrow=2, ncol=2),
                                       matrix(nrow=3, ncol=2))), 2:3)
})

test_that(".snippet", {
    L <- paste0(LETTERS[1:26], collapse="")
    l <- paste0(letters[1:26], collapse="")
    expect_equal(topdownr:::.snippet(L, 100), L)
    expect_equal(topdownr:::.snippet(L, 10), "ABCD...XYZ")
    expect_equal(topdownr:::.snippet(L, 11), "ABCD...WXYZ")
    expect_equal(topdownr:::.snippet(c(l, L), 10), c("abcd...xyz",
                                                    "ABCD...XYZ"))
    expect_equal(topdownr:::.snippet(c(l, L), 11), c("abcd...wxyz",
                                                    "ABCD...WXYZ"))
})

test_that(".subset", {
    expect_error(topdownr:::.subset(1:2, 10, letters[1:2]))
    expect_error(topdownr:::.subset(c(1, NA, 2), 10, letters[1:10]),
                 "'NA' is not supported")
    expect_error(topdownr:::.subset(list(foo=1:10), 10, letters[1:10]),
                 "Unknown")
    expect_equal(topdownr:::.subset(1:2, 10, letters[1:10]), 1:2)
    expect_equal(topdownr:::.subset(Rle(1:2, c(1, 1)), 10, letters[1:10]), 1:2)
    expect_equal(topdownr:::.subset(c(TRUE, TRUE, rep(FALSE, 8)), 10,
                                   letters[1:10]), 1:2)
    expect_equal(topdownr:::.subset(c("a", "b"), 10, letters[1:10]), 1:2)
})

test_that(".subsetByCharacter", {
    expect_error(topdownr:::.subsetByCharacter(1:2, LETTERS[1:2]))
    expect_error(topdownr:::.subsetByCharacter(letters[1:2], TRUE))
    expect_error(topdownr:::.subsetByCharacter(letters[1:2], LETTERS[1:2]),
                 "Subscript out of bound: 'a', 'b'")
    expect_equal(topdownr:::.subsetByCharacter(letters[1:2], letters[4:1]), 4:3)
    expect_equal(topdownr:::.subsetByCharacter(letters[1:2]), integer())
})

test_that(".subsetByLogical", {
    expect_error(topdownr:::.subsetByLogical(1:2, 10))
    expect_error(topdownr:::.subsetByLogical(TRUE, TRUE))
    expect_error(topdownr:::.subsetByLogical("foo", 10))
    expect_equal(topdownr:::.subsetByLogical(TRUE, 10), 1:10)
    expect_equal(topdownr:::.subsetByLogical(c(TRUE, FALSE), 10),
                 seq(1, 10, by=2))
    expect_equal(topdownr:::.subsetByLogical(rep(TRUE, 10), 10), 1:10)
    expect_equal(topdownr:::.subsetByLogical(rep(TRUE, 12), 10), 1:10)
})

test_that(".subsetByNumeric", {
    expect_error(topdownr:::.subsetByNumeric(TRUE, 10))
    expect_error(topdownr:::.subsetByNumeric(1:10, TRUE))
    expect_error(topdownr:::.subsetByNumeric("foo", 10))
    expect_error(topdownr:::.subsetByNumeric(c(1, 3, 12), 10),
                 "Subscript out of bound: '12'")
    expect_equal(topdownr:::.subsetByNumeric(1:10, 20), 1:10)
})

test_that(".subsetFiles", {
    expect_equal(topdownr:::.subsetFiles(
                  c("foo.experiments.csv", "foo.mzML", "bar.txt"), "foo"),
                 c(TRUE, TRUE, FALSE))
})

test_that(".swapFileExt", {
    expect_equal(topdownr:::.swapFileExt("foo.xml"), "foo.meth")
    expect_equal(topdownr:::.swapFileExt("foo.xml", "bar"), "foo.bar")
})

test_that(".targetedMassListToMz", {
    expect_error(topdownr:::.targetedMassListToMz(1:3))
    expect_equal(topdownr:::.targetedMassListToMz(c("(mz=1000.12 z=2 name=foo)",
                                                   "(mz=933.99 z=3 name=)")),
                 c(1000.1, 933.9))
})

test_that(".topDownFileExtRx", {
    ext <- c("experiments\\.csv", "fasta", "mz[Mm][Ll]", "raw", "txt")
    gz <- "(\\.(gz|bz2|xz|zip))?$"
    expect_error(topdownr:::.topDownFileExtRx("foo"))
    expect_equal(topdownr:::.topDownFileExtRx(),
                 paste0("\\.", ext[-4], gz, collapse="|"))
    expect_equal(topdownr:::.topDownFileExtRx("all"),
                 paste0("\\.", ext, gz, collapse="|"))
    expect_equal(topdownr:::.topDownFileExtRx("cfmt"),
                 paste0("\\.", ext[-4], gz, collapse="|"))
    expect_equal(topdownr:::.topDownFileExtRx("csv"),
                 paste0("\\.", ext[1], gz, collapse="|"))
    expect_equal(topdownr:::.topDownFileExtRx("mzml"),
                 paste0("\\.", ext[3], gz, collapse="|"))
    expect_equal(topdownr:::.topDownFileExtRx("txt"),
                 paste0("\\.", ext[5], gz, collapse="|"))
})

test_that(".topIdx", {
    d <- 1:10
    g <- rep_len(LETTERS[1:3], 10)
    expect_error(topdownr:::.topIdx(logical(10)),
                 "'x' has to be of type")
    expect_error(topdownr:::.topIdx(d, groupBy=g, n=-1),
                 "'n' has to be greater or equal than 1.")
    expect_error(topdownr:::.topIdx(d, groupBy=1:3, n=3), "have to be equal.")
    expect_equal(topdownr:::.topIdx(d, groupBy=g, n=3),
                 c(10, 7, 4, 8, 5, 2, 9, 6, 3))
    expect_equal(topdownr:::.topIdx(d, groupBy=g, n=2),
                 c(10, 7, 8, 5, 9, 6))
})

test_that(".translateThermoIdToScanId", {
    expect_error(topdownr:::.translateThermoIdToScanId(NULL))
    expect_error(topdownr:::.translateThermoIdToScanId(1:10))
    expect_error(topdownr:::.translateThermoIdToScanId(c("", "")))
    expect_equal(topdownr:::.translateThermoIdToScanId(c(
        "controllerType=0 controllerNumber=1 scan=11",
        "controllerType=0 controllerNumber=1 scan=12",
        "controllerType=0 controllerNumber=1 scan=13",
        "controllerType=0 controllerNumber=1 scan=14",
        "controllerType=0 controllerNumber=1 scan=15",
        "controllerType=0 controllerNumber=1 scan=16")), 11:16)
    expect_equal(topdownr:::.translateThermoIdToScanId(c(
        "scan=21 file=191",
        "scan=22 file=191",
        "scan=23 file=191",
        "scan=24 file=191",
        "scan=25 file=191",
        "scan=26 file=191")), 21:26)
})
