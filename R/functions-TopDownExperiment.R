#' TopDownExperiment
#'
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
TopDownExperiment <- function(path, pattern=".*",
                              type=c("a", "b", "c", "x", "y", "z"),
                              modifications=c(C=57.02146),
                              neutralLoss=defaultNeutralLoss(),
                              tolerance=10e-6,
                              verbose=interactive(), ...) {
  tdf <- .readTopDownFiles(path=path, pattern=pattern, verbose=verbose)

  sequence <- tdf$fasta

  header <- .mergeScanConditionAndHeaderInformation(tdf$ScanConditions,
                                                    tdf$HeaderInformation)
  msnexp <- .mergeSpectraAndHeaderInformation(tdf$MSnExp, header)

  ftab <- .calculateFragments(sequence, type=type, modifications=modifications,
                              neutralLoss=neutralLoss, verbose=verbose)

  atab <- .assignmentTable(msnexp, ftab, tolerance=tolerance,
                           verbose=verbose, ...)

  td <- new("TopDownExperiment",
            assayData=assayData(msnexp),
            featureData=featureData(msnexp),
            phenoData=phenoData(msnexp),
            experimentData=experimentData(msnexp),
            processingData=processingData(msnexp),
            sequence=sequence,
            fragmentTable=ftab,
            assignmentTable=atab,
            ...)

  if (validObject(td)) {
    td
  }
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
  n <- .vapply1d(a, function(aa)length(aa[[1L]]))
  data.table(SpectrumId=rep.int(featureNames(msnexp), n),
             FragmentId=as.double(unlist(lapply(a, "[[", "FragId"))),
             MzId=as.double(unlist(lapply(a, "[[", "MzId"))),
             key=c("SpectrumId", "FragmentId", "MzId"))
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

#' Get fragment types.
#'
#' @param object TopDownExperiment
#' @return character
#' @noRd
.fragmentTypes <- function(object) {
  fragmentTable(object)$type
}
