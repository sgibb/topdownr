#' Subset MSnbase::Spectrum2 objects
#'
#' @param object Spectrum2 object
#' @param i double, id of mz values
#' @return subset of Spectrum2
#' @noRd
.subsetSpectrum2 <- function(object, i) {
  object@mz <- mz(object)[i]
  object@intensity <- intensity(object)[i]
  object@tic <- sum(intensity(object))
  object@peaksCount <- length(mz(object))
  if (validObject(object)) {
    return(object)
  }
}
