#' TopDownExperiment
#'
#' @param sequence \code{character}, amino acid sequence.
#' @param path \code{character}, path to top down files.
#' @param pattern \code{character}, pattern for file names.
#' @param type character, type of fragments (should be some more than in
#' MSnbase)
#' @param modification modifications (see MSnbase::calculateFragments)
#' @param neutralLoss neutral loss (see MSnbase::calculateFragments)
#' @param tolerance double, tolerance to match peaks
#' @param verbose logical, verbose output?
#' @return TopDownExperiment object.
#' @export
TopDownExperiment <- function(sequence, path, pattern,
                              type=c("a", "b", "c", "x", "y", "z"),
                              modifications=c(C=57.02146),
                              neutralLoss=defaultNeutralLoss(),
                              tolerance=10e-6,
                              verbose=interactive(), ...) {
  tdf <- .readTopDownFiles(path=path, pattern=pattern, verbose=verbose)
  header <- .mergeScanConditionAndHeaderInformation(tdf$ScanConditions,
                                                    tdf$HeaderInformation)
  msnexp <- .mergeSpectraAndHeaderInformation(tdf$MSnExp, header)

  ftab <- .calculateFragments(sequence, type=type, modifications=modifications,
                              neutralLoss=neutralLoss, verbose=verbose)

  atab <- .assignmentTable(msnexp, ftab, tolerance=tolerance,
                           verbose=verbose, ...)

  new("TopDownExperiment",
      sequence=sequence,
      msnExp=msnexp,
      fragmentTable=ftab,
      assignmentTable=atab,
      ...)
}

#' Create assignment table/match peaks.
#'
#' @param msnexp MSnExp object
#' @param ftab data.table, theoreticalFragmentTable fragments
#' @param tolerance double, tolerance to match peaks
#' @param verbose logical, verbose output?
#' @param \ldots further arguments passed to MSnbase::spectrapply
#' @return data.table
#' @noRd
.assignmentTable <- function(msnexp, ftab, tolerance=25e-6,
                             verbose=interactive(), ...) {
  .msg(verbose, "Looking for fragments in spectra")
  a <- spectrapply(msnexp, function(s) {
    i <- MSnbase:::matchPeaks(s, y=ftab$mz, tolerance=tolerance)
    notNA <- !is.na(i)
    list(FragId=i[notNA], MzId=which(notNA))
  }, ...)
  n <- vapply(a, function(aa)length(aa[[1L]]), double(1L))
  data.table(SpectrumId=rep.int(as.double(fData(msnexp)$spectrum), n),
             FragmentId=as.double(unlist(lapply(a, "[[", "FragId"))),
             MzId=as.double(unlist(lapply(a, "[[", "MzId"))),
             key=c("SpectrumId", "FragmentId", "MzId"))
}

#' Filter by fragment type.
#'
#' @param object TopDownExperiment object
#' @param type character, fragment type
#' @return subsetted object
#' @noRd
.filterFragmentType <- function(object, type) {
  # seems to be faster than ftab[,unique(type)]
  types <- .fragmentTypes(object)
  utypes <- unique(types)
  if (!all(type %in% utypes)) {
    stop("Type ",
         paste0(dQuote(type[!type %in% utypes]), collapse=", "),
         " is not valid!")
  }
  fid <- .fragmentId(object)[types %in% type]
  atab <- assignmentTable(object)[FragmentId %in% fid, ]

  object@msnExp <- .subsetMSnExpSpectra(msnExp(object), atab$SpectrumId,
                                        atab$MzId)
  object@assignmentTable <- atab

  if (validObject(object)) {
    object
  }
}

#' Get fragment ID.
#'
#' @param object
#' @return numeric
.fragmentId <- function(object) {
  fragmentTable(object)$FragmentId
}

#' Get fragment types.
#'
#' @param object
#' @return character
.fragmentTypes <- function(object) {
  fragmentTable(object)$type
}
