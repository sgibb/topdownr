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
      assayData=msnexp@assayData,
      phenoData=msnexp@phenoData,
      featureData=msnexp@featureData,
      processingData=msnexp@processingData,
      experimentData=msnexp@experimentData,
      .cache=msnexp@.cache,
      fragmentTable=ftab,
      assignmentTable=atab,
      ...)
}

#' create assignment table / match peaks
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
    i[!is.na(i)]
  }, ...)
  data.table(SpectrumId=rep.int(as.double(fData(msnexp)$spectrum), lengths(a)),
             FragmentId=as.double(unlist(a)),
             key=c("SpectrumId", "FragmentId"))
}
