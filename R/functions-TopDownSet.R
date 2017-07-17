#' Read TopDown files.
#'
#' Read all TopDown files:
#'  - .fasta (peptide sequence)
#'  - .mzML (spectra)
#'  - .experiments.csv (fragmentation conditions)
#'  - .txt (header information)
#'
#' This is the regular \code{\linkS4class{TopDownSet}} constructor.
#'
#' @param path character, file path
#' @param pattern character, filename pattern
#' @param onDisk logical, return MSnExp or (if TRUE) OnDiskMSnExp
#' @param verbose logical, verbose output?
#' @return list (splitted by file extension) with file path
#' @export
#' @noRd
readTopDownFiles <- function(path, pattern=".*",
                             type=c("a", "b", "c", "x", "y", "z"),
                             modifications=c(C=57.02146),
                             neutralLoss=defaultNeutralLoss(),
                             tolerance=10e-6,
                             verbose=interactive(), ...) {

  files <- .listTopDownFiles(path, pattern=pattern)

  if (any(lengths(files)) == 0L) {
    ext <- c("experiments.csv", "fasta", "mzML", "txt")
    stop("Could not found any ", paste0(ext[lengths(files) == 0L],
                                        collapse=" or "), " files!")
  }

  sequence <- .readFasta(files$fasta, verbose=verbose)

  fragmentViews <- .calculateFragments(sequence=sequence,
                                       type=type,
                                       modifications=modifications,
                                       neutralLoss=neutralLoss)

  scanConditions <- do.call(rbind, lapply(files$txt, .readScanHeadsTable,
                                          verbose=verbose))

  headerInformation <- do.call(rbind, lapply(files$csv, .readExperimentCsv,
                                             verbose=verbose))

  mzml <- mapply(.readMzMl,
                 file=files$mzML,
                 scans=split(scanConditions$Scan, scanConditions$File),
                 MoreArgs=list(fmass=elementMetadata(fragmentViews)$mass,
                               tolerance=tolerance),
                 SIMPLIFY=FALSE)

  mzmlHeader <- do.call(rbind, lapply(mzml, "[[", "hd"))

  scanHeadsman <- .mergeScanConditionAndHeaderInformation(scanConditions,
                                                          headerInformation)

  header <- .mergeSpectraAndHeaderInformation(mzmlHeader, scanHeadsman)

  assay <- do.call(cbind, lapply(mzml, "[[", "m"))
  dimnames(assay) <- list(names(fragmentViews),
                          rownames(header))

  new("TopDownSet",
      rowViews=fragmentViews,
      colData=.colsToRle(as(header, "DataFrame")),
      assay=assay,
      files=basename(unlist(unname(files))),
      processing=.logmsg("Data loaded."))
}

#' Test for TopDownSet class
#'
#' @param object object to test
#' @return TRUE if object is a TopDownSet otherwise fails with an error
#' @noRd
.isTopDownSet <- function(object) {
  if (!isTRUE(is(object, "TopDownSet"))) {
    stop("'object' has to be an 'TopDownSet' object.")
  }
  TRUE
}

#' @noRd
fragmentMass <- function(object) {
  .isTopDownSet(object)
  elementMetadata(object@rowViews)$mass
}

#' @noRd
fragmentNames <- function(object) {
  .isTopDownSet(object)
  names(object@rowViews)
}

#' @noRd
fragmentType <- function(object) {
  .isTopDownSet(object)
  elementMetadata(object@rowViews)$type
}

#' Create NCB Map (N-/C-terminal, or both)
#'
#' @param object TopDownSet
#' @return Matrix, Nterm == 1, Cterm == 2, both == 3
#' @noRd
.ncbMap <- function(object, nterm=c("a", "b", "c"), cterm=c("x", "y", "z")) {
  .isTopDownSet(object)

  w <- width(object@rowViews)
  mn <- mc <- object@assay
  selN <- fragmentType(object) %in% nterm
  selC <- fragmentType(object) %in% cterm
  mn[!selN,] <- 0L
  mc[!selC,] <- 0L

  mn <- as(.colSumsGroup(mn, w) > 0L, "dgCMatrix")
  mc <- as(.colSumsGroup(mc, max(w) + 1L - w) > 0L, "dgCMatrix")
  mc@x[] <- 2
  mn + mc
}

#' @noRd
.tdsLogMsg <- function(object, msg) {
  .isTopDownSet(object)
  object@processing <- c(object@processing, .logmsg(msg))
  object
}

#' Validate TopDownSet
#'
#' @param object TopDownSet
#' @return TRUE (if valid) else character with msg what was incorrect
#' @noRd
.validateTopDownSet <- function(object) {
  msg <- character()

  if (nrow(object@assay) != length(object@rowViews)) {
    msg <- c(msg, "Mismatch between fragment data in 'rowViews' and 'assay'.")
  }

  if (any(rownames(object@assay) != names(object@rowViews))) {
    msg <- c(msg,
             "Mismatch between fragment names in 'rowViews' and 'assay'.")
  }

  if (ncol(object@assay) != nrow(object@colData)) {
    msg <- c(msg, "Mismatch between condition data in 'colData' and 'assay'.")
  }

  if (any(colnames(object@assay) != rownames(object@colData))) {
    msg <- c(msg,
             "Mismatch between condition names in 'colData' and 'assay'.")
  }

  if (length(msg)) {
    msg
  } else {
    TRUE
  }
}
