#' TopDownExperiment
#'
#' @param path \code{character}, path to top down files.
#' @param pattern \code{character}, pattern for file names.
#' @param type character, type of fragments (should be some more than in
#' MSnbase)
#' @param modifications modifications (see MSnbase::calculateFragments)
#' @param neutralLoss neutral loss (see MSnbase::calculateFragments)
#' @param tolerance double, tolerance to match peaks
#' @param onDisk logical, use \code{\linkS4class{MSnExp}} (\code{FALSE},
#' default) or \code{\linkS4class{OnDiskMSnExp}} (\code{TRUE}) for spectra
#' storage.
#' @param verbose logical, verbose output?
#' @return TopDownExperiment object.
#' @export
TopDownExperiment <- function(path, pattern=".*",
                              type=c("a", "b", "c", "x", "y", "z"),
                              modifications=c(C=57.02146),
                              neutralLoss=defaultNeutralLoss(),
                              tolerance=10e-6,
                              onDisk=FALSE,
                              verbose=interactive(), ...) {
  tdf <- .readTopDownFiles(path=path, pattern=pattern, onDisk=onDisk,
                           verbose=verbose)

  sequence <- tdf$fasta

  header <- .mergeScanConditionAndHeaderInformation(tdf$ScanConditions,
                                                    tdf$HeaderInformation)
  msnexp <- .mergeSpectraAndHeaderInformation(tdf$MSnExp, header)
  msnexp@processingData@processing <- character()

  msnexp <- .logmsg(msnexp, "Data loaded.")

  ftab <- .calculateFragments(sequence, type=type, modifications=modifications,
                              neutralLoss=neutralLoss, verbose=verbose)

  m <- .matchFragments(msnexp, ftab, tolerance=tolerance, verbose=verbose, ...)

  td <- new("TopDownExperiment",
            assayData=m$assayData,
            featureData=featureData(msnexp),
            phenoData=phenoData(msnexp),
            experimentData=experimentData(msnexp),
            processingData=processingData(msnexp),
            sequence=sequence,
            fragmentTable=ftab,
            assignmentTable=m$assignmentTable,
            ...)

  pc <- c(sum(peaksCount(td)), sum(peaksCount(msnexp)))
  msg <- sprintf("Fragments matched: %d/%d (%d %%) peaks kept.",
                 pc[1L], pc[2L], round(pc[1L]/pc[2L] * 100L))
  td <- .logmsg(td, msg)
  .msg(verbose, msg)

  if (validObject(td)) {
    td
  }
}

#' Filter by fragment id.
#'
#' @param object TopDownExperiment object
#' @param id double, fragment id
#' @return subsetted object
#' @noRd
.filterFragmentId <- function(object, id) {
  atab <- assignmentTable(object)[FragmentId %in% id, ]
  if (!nrow(atab)) {
    stop("Selecting would result in an empty object. No subsetting applied.")
  }
  object <- .subsetMSnExpSpectra(object, atab$SpectrumId, atab$MzId)
  object@assignmentTable <- .updateAssignmentTableMzId(atab)

  if (validObject(object)) {
    object
  }
}

#' Filter by fragment ion or type.
#'
#' @param object TopDownExperiment object
#' @param ion character, fragment ion
#' @return subsetted object
#' @noRd
.filterFragmentIonOrType <- function(object, ion) {
  # seems to be faster than ftab[,unique(ion)]
  ions <- .fragmentIons(object)
  types <- .fragmentTypes(object)
  u <- unique(c(ions, types))
  if (!all(ion %in% u)) {
    stop("Ion(s)/Type(s) ",
         paste0(dQuote(ion[!ion %in% u]), collapse=", "),
         " not found!")
  }
  .filterFragmentId(object, .fragmentId(object)[ions %in% ion | types %in% ion])
}

#' Filter by fragment pos.
#'
#' @param object TopDownExperiment object
#' @param pos integer, fragment/bound position
#' @return subsetted object
#' @noRd
.filterFragmentPos <- function(object, pos) {
  # seems to be faster than ftab[,unique(pos)]
  positions <- .fragmentPos(object)
  uPos <- unique(positions)
  if (!all(pos %in% uPos)) {
    stop("Position(s) ",
         paste0(dQuote(pos[!pos %in% uPos]), collapse=", "),
         " not found!")
  }
  .filterFragmentId(object, .fragmentId(object)[positions %in% pos])
}

#' Get fragment ID.
#'
#' @param object TopDownExperiment
#' @return numeric
#' @noRd
.fragmentId <- function(object) {
  fragmentTable(object)$FragmentId
}

#' Get fragment ions.
#'
#' @param object TopDownExperiment
#' @return character
#' @noRd
.fragmentIons <- function(object) {
  fragmentTable(object)$ion
}

#' Get fragment positions.
#'
#' @param object TopDownExperiment
#' @return integer
#' @noRd
.fragmentPos <- function(object) {
  fragmentTable(object)$pos
}

#' Get fragment types.
#'
#' @param object TopDownExperiment
#' @return character
#' @noRd
.fragmentTypes <- function(object) {
  fragmentTable(object)$type
}

#' Add log message.
#'
#' @param object TopDownExperiment
#' @param msg character, log message
#' @return TopDownExperiment
#' @noRd
.logmsg <- function(object, msg, date=TRUE) {
  if (date) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  }
  object@processingData@processing <-
    c(object@processingData@processing, msg)
  object
}

#' Look for fragment/mz matching.
#'
#' Removes non-matching peaks and creates assignment table.
#'
#' @param msnexp MSnExp object
#' @param ftab data.table, theoreticalFragmentTable fragments
#' @param tolerance double, tolerance to match peaks
#' @param verbose logical, verbose output?
#' @param \ldots further arguments passed to MSnbase::spectrapply
#' @return list, new assayData (environment) and assignmentTable (data.table)
#' @noRd
.matchFragments <- function(msnexp, ftab, tolerance=25e-6,
                            verbose=interactive(), ...) {
  .msg(verbose, "Looking for fragments in spectra.")

  if (verbose) {
    pb <- txtProgressBar(min=0L, max=length(msnexp), style=3L)
    on.exit(close(pb))
    i <- 0L
  }

  newAssay <- new.env(parent=emptyenv())
  fragId <- vector(mode="list", length=length(msnexp))
  names(fragId) <- featureNames(msnexp)
  oldAssay <- assayData(msnexp)

  for (specName in featureNames(msnexp)) {
    if (verbose) {
      i <- i + 1L
      setTxtProgressBar(pb, i)
    }
    sp <- get(specName, oldAssay)
    m <- MSnbase:::matchPeaks(sp, y=ftab$mz, tolerance=tolerance)
    notNA <- which(!is.na(m))
    sp <- .subsetSpectrum2(sp, notNA)
    fragId[[specName]] <- m[notNA]
    assign(specName, sp, newAssay)
  }

  lockEnvironment(newAssay, bindings=TRUE)

  d <- data.table(SpectrumId=rep.int(featureNames(msnexp), lengths(fragId)),
                  FragmentId=as.double(unlist(fragId)),
                  key=c("SpectrumId", "FragmentId"))
  d <- .updateAssignmentTableMzId(d)

  list(assayData=newAssay, assignmentTable=d)
}

#' Validate TopDownExperiment
#'
#' @param object TopDownExperiment
#' @return TRUE (if valid) else character with msg what was incorrect
#' @noRd
.validateTopDownExperiment <- function(object) {
  msg <- character()

  if (!all(object@assignmentTable$SpectrumId %in% featureNames(object))) {
    msg <- c(msg, "IDs in assignment table don't match feature names.")
  }

  if (any(unlist(lapply(peaksCount(object), seq_len), use.names=FALSE) !=
          object@assignmentTable$MzId)) {
    msg <- c(msg, "Mismatch in spectra and assignment table's peak indices.")
  }

  if (!all(object@assignmentTable$FragmentId %in% object@fragmentTable$FragmentId)) {
    msg <- c(msg, "Mismatch in fragment and assignment table's fragment indices.")
  }

  if (is.unsorted(object@fragmentTable$mz)) {
    msg <- c(msg, "Mz values have to be sorted in the fragment table.")
  }

  if (length(msg)) {
    msg
  } else {
    TRUE
  }
}
