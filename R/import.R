#' Test for file existence with error message.
#' @param file filename
#' @noRd
.fileExists <- function(file) {
  if (!file.exists(file)) {
    stop(file, " doesn't exists!")
  }
  TRUE
}

#' List TopDown files
#'
#' List all TopDown files:
#'  - .mzML (spectra)
#'  - .experiments.csv (fragmentation conditions)
#'  - .txt (header information)
#'
#' @param path file path
#' @param pattern filename pattern
#' @return list (splitted by file extension) with file path
#' @noRd
.listTopDownFiles <- function(path, pattern=".*") {
  files <- list.files(path,
                      pattern=paste0("(", pattern, ")",
                                     c("\\.mzML$", "\\.experiments\\.csv$",
                                       "\\.txt$"), collapse="|"),
                      full.names=TRUE)
  split(files, file_ext(files))
}

#' Read TopDown files.
#'
#' Read all TopDown files:
#'  - .mzML (spectra)
#'  - .experiments.csv (fragmentation conditions)
#'  - .txt (header information)
#'
#' @param path character, file path
#' @param pattern character, filename pattern
#' @param verbose logical, verbose output?
#' @return list (splitted by file extension) with file path
#' @noRd
.readTopDownFiles <- function(path, pattern=".*", verbose=interactive()) {
  files <- .listTopDownFiles(path, pattern=pattern)

  if (any(lengths(files)) == 0L) {
    stop("Could not found any ",
         paste0(c("experiments.csv", "mzML", "txt")[lengths(files) == 0L],
                collapse=" or "), " files!")
  }

  mzml <- .readMSData2(files$mzML, verbose=verbose)
  csv <- do.call(rbind, lapply(files$csv, .readExperimentCsv,
                               verbose=verbose))
  txt <- do.call(rbind, lapply(files$txt, .readScanHeadsTable,
                               verbose=verbose))
  list(MSnExp=mzml, ScanConditions=csv, HeaderInformation=txt)
}

#' Read ScanHeadMans method (experiments.csv) output.
#'
#' 1 row per condition
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return data.table
#' @noRd
.readExperimentCsv <- function(file, verbose=interactive()) {
  stopifnot(file_ext(file) == "csv")
  d <- fread(file, showProgress=verbose)

  .msg(verbose, "Reading ", nrow(d), " experiment conditions from file ",
       basename(file))

  ## drop MS1
  d <- d[MSLevel == 2L,]

  d[, ConditionId := .I]
  d[, Mz := .targetedMassListToMz(TargetedMassList)]
  d[, File := basename(file)]
}

#' Read ScanHeadMans header (txt) output.
#'
#' 1 row per scan (could have more rows than experiment.csv, because multiple
#' scans per condition are possible).
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return data.table
#' @noRd
.readScanHeadsTable <- function(file, verbose=interactive()) {
  stopifnot(file_ext(file) == "txt")
  d <- fread(file, showProgress=verbose)

  .msg(verbose, "Reading ", nrow(d), " header information from file ",
       basename(file))

  ## drop MS1
  d <- d[MSOrder == 2L,]

  d[, File := basename(file)]
}

#' Read MS2 Spectra (mzML)
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return MSnExp
#' @noRd
.readMSData2 <- function(files, verbose=interactive()) {
  readMSData2(files, msLevel.=2, verbose=verbose)
}