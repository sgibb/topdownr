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
                      recursive=TRUE, full.names=TRUE)
  l <- split(files, .fileExt(files))
  n <- lengths(l)

  ext <- c("csv", "fasta", "mzML", "txt")

  if (!length(n) || any(!ext %in% names(l))) {
    stop("Could not find any ", paste0(ext[!ext %in% names(l)],
                                       collapse=", "), " files!")
  }

  if (n["fasta"] > 1L) {
    stop("More than one fasta file found. Consider the 'pattern' argument.")
  }

  if (!all(n["csv"] == n[!grepl("fasta", names(n))])) {
    nd <- n[!grepl("fasta", names(n))]
    stop("There have to be the same number of csv, mzML and txt files. ",
         "Found: ", paste0(names(nd), "=", nd, collapse=", "))
  }
  l
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
  stopifnot(.fileExt(file) == "csv")
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
  stopifnot(.fileExt(file) == "txt")
  d <- read.csv(file, stringsAsFactors=FALSE)
  colnames(d) <- .formatNames(colnames(d))

  .msg(verbose, "Reading ", nrow(d), " header information from file ",
       basename(file))

  ## drop MS1
  d <- d[d$MSOrder == 2L,]

  # TODO: somehow the FilterString doesn't always contains the right mass label.
  # For now we just take the first non-duplicated (unique) condition
  # disabled: d$Condition <- as.integer(.filterStringToId(d$FilterString))
  #
  # See the following issues for details:
  # - https://github.com/sgibb/topdown/issues/14
  # - https://github.com/sgibb/topdown/issues/25
  d$FilterString <- .fixFilterStringId(d$FilterString)
  d$Condition <- cumsum(!duplicated(d$FilterString))

  d[is.na(d)] <- 0L

  d$ETDActivation[d$Activation1 == "ETD"] <- d$Energy1[d$Activation1 == "ETD"]
  d$CIDActivation[d$Activation1 == "CID"] <- d$Energy1[d$Activation1 == "CID"]
  d$CIDActivation[d$Activation2 == "CID"] <- d$Energy2[d$Activation2 == "CID"]
  d$HCDActivation[d$Activation1 == "HCD"] <- d$Energy1[d$Activation1 == "HCD"]
  d$HCDActivation[d$Activation2 == "HCD"] <- d$Energy2[d$Activation2 == "HCD"]

  d[is.na(d)] <- 0L

  d$Activation <- .fragmentationMethod(d[, paste0(c("ETD", "CID", "HCD"),
                                                  "Activation")])

  d$ActivationString <- paste(.formatNumbers(d$ETDActivation),
                              .formatNumbers(d$CIDActivation),
                              .formatNumbers(d$HCDActivation), sep=":")

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
  hd <- hd[i, !colnames(hd) %in% c("injectionTime", "seqNum"), drop=FALSE]
  colnames(hd)[grepl("acquisitionNum", colnames(hd), fixed=TRUE)] <- "Scan"
  hd$File <- gsub(.topDownFileExtRx("mzml"), "", basename(file))
  colnames(hd) <- .formatNames(colnames(hd))

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
                        round(sum(m != 0L)/sum(hd$PeaksCount) * 100, 1L)))

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
  merge(mzml, scdm, sort=FALSE, by=c("File", "Scan"),
        suffixes=c(".SpectraInformation", ".HeaderInformation"))
}
