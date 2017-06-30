#' Subset MSnExp by spectra and the latter by mz
#'
#' @param msnexp MSnExp object
#' @param spectrumId double, id of spectra
#' @param mzId double, ids of mz values
#' @return subset of MSnExp
#' @noRd
.subsetMSnExpSpectra <- function(msnexp, spectrumId, mzId) {
  stopifnot(length(spectrumId) == length(mzId))
  msnexp <- msnexp[unique(spectrumId)]
  msnexp@assayData <- as.environment(
    mapply(.subsetSpectrum2,
           object=mget(featureNames(msnexp), assayData(msnexp)),
           i=split(mzId, spectrumId), SIMPLIFY=FALSE))
  if (validObject(msnexp)) {
    msnexp
  }
}
