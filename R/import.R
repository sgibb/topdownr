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
#'  - .fasta (peptide sequence)
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
                      pattern=paste0(pattern, "(",
                                     .topDownFileExtRx("cfmt"), ")"),
                      full.names=TRUE)
  split(files, file_ext(files))
}

#' Read TopDown files.
#'
#' Read all TopDown files:
#'  - .fasta (peptide sequence)
#'  - .mzML (spectra)
#'  - .experiments.csv (fragmentation conditions)
#'  - .txt (header information)
#'
#' @param path character, file path
#' @param pattern character, filename pattern
#' @param onDisk logical, return MSnExp or (if TRUE) OnDiskMSnExp
#' @param verbose logical, verbose output?
#' @return list (splitted by file extension) with file path
#' @noRd
.readTopDownFiles <- function(path, pattern=".*", onDisk=FALSE,
                              verbose=interactive()) {
  files <- .listTopDownFiles(path, pattern=pattern)

  if (any(lengths(files)) == 0L) {
    ext <- c("experiments.csv", "fasta", "mzML", "txt")
    stop("Could not found any ", paste0(ext[lengths(files) == 0L],
                                        collapse=" or "), " files!")
  }

  fasta <- .readFasta(files$fasta, verbose=verbose)
  mzml <- .readMSData(files$mzML, onDisk=onDisk, verbose=verbose)
  csv <- do.call(rbind, lapply(files$csv, .readExperimentCsv,
                               verbose=verbose))
  txt <- do.call(rbind, lapply(files$txt, .readScanHeadsTable,
                               verbose=verbose))
  list(fasta=fasta, MSnExp=mzml, ScanConditions=csv, HeaderInformation=txt)
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
  d <- read.csv(file, stringsAsFactors=FALSE)
  colnames(d) <- .formatNames(colnames(d))

  .msg(verbose, "Read ", nrow(d), " experiment conditions from file ",
       basename(file))

  ## drop MS1
  d <- d[d$MSLevel == 2L,]

  d[is.na(d)] <- 0L

  d$ConditionId <- seq_len(nrow(d))
  d$Mz <- .targetedMassListToMz(d$TargetedMassList)
  d$File <- gsub(.topDownFileExtRx("csv"), "", basename(file))
  d
}

#' Read fasta file.
#'
#' Stupid fasta reader, ignores comments, and keeps just the first sequence.
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return character
#' @noRd
.readFasta <- function(file, verbose=interactive()) {
  aa <- readAAStringSet(file, nrec=1L, use.names=FALSE)[[1L]]
  .msg(verbose, "Read sequence from fasta file ", basename(file))
  if (!length(aa)) {
    stop("No sequence found.")
  }
  aa
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
  colnames(d) <- .formatNames(colnames(d))

  .msg(verbose, "Reading ", nrow(d), " header information from file ",
       basename(file))

  ## drop MS1
  d <- d[MSOrder == 2L,]

  d[Activation1 == "ETD", ETDActivation := Energy1]
  d[Activation1 == "CID", CIDActivation := Energy1]
  d[Activation2 == "CID", CIDActivation := Energy2]
  d[Activation1 == "HCD", HCDActivation := Energy1]
  d[Activation2 == "HCD", HCDActivation := Energy2]

  d[is.na(d)] <- 0L

  d[, ConditionId := as.integer(.filterStringToId(FilterString))]

  d[, File := gsub(.topDownFileExtRx("txt"), "", basename(file))]
}

#' Read MS2 Spectra (mzML)
#'
#' @param file character, filename
#' @param onDisk logical, return MSnExp or (if TRUE) OnDiskMSnExp
#' @param verbose logical, verbose output?
#' @return MSnExp
#' @noRd
.readMSData <- function(files, onDisk=FALSE, verbose=interactive()) {
  if (onDisk) {
    readMSData2(files, msLevel.=2L, verbose=verbose)
  } else {
    readMSData(files, msLevel.=2L, verbose=verbose)
  }
}

#' Merge ScanCondition and HeaderInformation
#'
#' @param sc data.table, scan conditions
#' @param hi data.table, header information
#' @return data.table
#' @noRd
.mergeScanConditionAndHeaderInformation <- function(sc, hi) {
  stopifnot(is(sc, "data.frame"))
  stopifnot(is(hi, "data.frame"))
  merge(sc, hi, by=c("File", "ConditionId"), all.y=TRUE,
        suffixes=c(".ScanCondition", ".HeaderInformation"))
}

#' Merge spectra and ScanConditions/HeaderInformation (into featureData slot)
#' @param msnexp MSnExp
#' @param header data.table, header information
#' @return modified MSnExp
#' @noRd
.mergeSpectraAndHeaderInformation <- function(msnexp, hi) {
  fd <- fData(msnexp)
  fd$File <- gsub(.topDownFileExtRx("mzml"), "",
                  basename(fileNames(msnexp))[fromFile(msnexp)])
  fd$Scan <- acquisitionNum(msnexp)
  m <- merge(fd, as.data.frame(hi), sort=FALSE)
  ## row.names are needed for fData<-
  sel <- match(m$spectrum, fd$spectrum)
  if (!length(sel)) {
    stop("The spectra and the header information have nothing in common.")
  }
  msnexp <- msnexp[sel]
  rownames(m) <- rownames(fd)[sel]
  fData(msnexp) <- m[, !grepl("^File$", colnames(m))]

  if(validObject(msnexp)) {
    return(msnexp)
  }
}
