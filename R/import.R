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

#' Read ScanHeadMans method (experiments.csv) output.
#'
#' 1 row per condition
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return data.frame
#' @noRd
.readExperimentCsv <- function(file, verbose=interactive()) {
  stopifnot(file_ext(file) == "csv")
  d <- read.csv(file, stringsAsFactors=FALSE)
  colnames(d) <- .formatNames(colnames(d))

  .msg(verbose, "Reading ", nrow(d), " experiment conditions from file ",
       basename(file))

  ## drop MS1
  d <- d[d$MSLevel == 2L,]

  d[is.na(d)] <- 0L

  d$Condition <- seq_len(nrow(d))
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
  .msg(verbose, "Reading sequence from fasta file ", basename(file))
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
#' @return data.frame
#' @noRd
.readScanHeadsTable <- function(file, verbose=interactive()) {
  stopifnot(file_ext(file) == "txt")
  d <- read.csv(file, stringsAsFactors=FALSE)
  colnames(d) <- .formatNames(colnames(d))

  .msg(verbose, "Reading ", nrow(d), " header information from file ",
       basename(file))

  ## drop MS1
  d <- d[d$MSOrder == 2L,]

  # TODO: somehow the FilterString doesn't always contains the right mass label.
  # For now we just take the first non-duplicated (unique) condition
  # disabled: d$Condition <- as.integer(.filterStringToId(d$FilterString))
  d$Condition <- cumsum(!duplicated(d$FilterString))

  d[is.na(d)] <- 0L

  d$ETDActivation[d$Activation1 == "ETD"] <- d$Energy1[d$Activation1 == "ETD"]
  d$CIDActivation[d$Activation1 == "CID"] <- d$Energy1[d$Activation1 == "CID"]
  d$CIDActivation[d$Activation2 == "CID"] <- d$Energy2[d$Activation2 == "CID"]
  d$HCDActivation[d$Activation1 == "HCD"] <- d$Energy1[d$Activation1 == "HCD"]
  d$HCDActivation[d$Activation2 == "HCD"] <- d$Energy2[d$Activation2 == "HCD"]

  d[is.na(d)] <- 0L

  d$File <- gsub(.topDownFileExtRx("txt"), "", basename(file))
  d
}

#' Read MS2 Spectra (mzML)
#'
#' @param file character, filename
#' @param fmass double, fragment mass
#' @param \ldots further arguments passed to .matchFragments
#' @param verbose logical, verbose output?
#' @return list (with data.frame for header and sparseMatrix with intensity
#' values)
#' @noRd
.readMzMl <- function(file, scans, fmass, ..., verbose=interactive()) {
  .msg(verbose, "Reading spectra information from file ", basename(file),
       appendLF=FALSE)

  fh <- openMSfile(file)
  on.exit(close(fh))

  hd <- header(fh)
  i <- which(hd$msLevel == 2L & hd$acquisitionNum %in% scans)
  hd <- hd[i, !grepl("seqNum", colnames(hd), fixed=TRUE), drop=FALSE]
  colnames(hd)[grepl("acquisitionNum", colnames(hd), fixed=TRUE)] <- "Scan"
  hd$File <- gsub(.topDownFileExtRx("mzml"), "", basename(file))

  nr <- nrow(hd)
  m <- Matrix(0L, nrow=length(fmass), ncol=nr, sparse=TRUE)

  for (j in seq_along(i)) {
    k <- .matchFragments(peaks(fh, i[j])[, 1L], fmass, ...)
    notNA <- !is.na(k)
    if (sum(notNA)) {
      m[k[notNA], j] <- peaks(fh, i[j])[notNA, 2L]
    }
  }

  .msg(verbose, sprintf(" (%02.1f%%)",
                        round(sum(m != 0L)/sum(hd$peaksCount) * 100, 1L)))

  list(hd=hd, m=m)
}

#' Merge ScanCondition and HeaderInformation
#'
#' @param sc data.frame, scan conditions
#' @param hi data.frame, header information
#' @return data.frame
#' @noRd
.mergeScanConditionAndHeaderInformation <- function(sc, hi) {
  stopifnot(is(sc, "data.frame"))
  stopifnot(is(hi, "data.frame"))
  merge(sc, hi, by=c("File", "Condition"), all.y=TRUE,
        suffixes=c(".ScanCondition", ".HeaderInformation"))
}

#' Merge spectra and ScanConditions/HeaderInformation (into featureData slot)
#'
#' @param mzml data.frame, header from mzML files
#' @param scdm data.frame, header from ScanHeadsman
#' @return merged data.frame
#' @noRd
.mergeSpectraAndHeaderInformation <- function(mzml, scdm) {
  merge(mzml, scdm, sort=FALSE)
}
