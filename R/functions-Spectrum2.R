#' aggregate MSnbase::Spectrum2 objects
#'
#' @param mz list, of mz values
#' @param int list, of intensity values
#' @param rt double, of rtime values
#' @param atab data.table, assignment table
#' @param ftab data.table, fragment table
#' @param \ldots arguments passed to new
#' @return subset of Spectrum2
#' @noRd
.aggregateSpectra <- function(mz, int, rt, atab, ftab, ...) {
  stopifnot(all(names(mz) == names(int)))
  stopifnot(all(atab$SpectrumId %in% names(mz)))
  stopifnot(all(lengths(mz) == lengths(int)))
  stopifnot(sum(lengths(mz)) == nrow(atab))
  m <- i <- matrix(NA_real_, nrow=nrow(ftab), ncol=length(mz),
                   dimnames=list(ftab$ions, names(mz)))
  coord <- cbind(atab$FragmentId, match(atab$SpectrumId, names(mz)))
  m[coord] <- unlist(mz)
  i[coord] <- unlist(int)
  mz <- rowMeans(m, na.rm=TRUE)
  int <- rowMeans(i, na.rm=TRUE)
  rt <- mean(rt, na.rm=TRUE)
  nan <- is.nan(mz)

  s <- new("Spectrum2", mz=mz[!nan], intensity=int[!nan], rt=rt, ...)

  if (validObject(s)) {
    return(s)
  }
}

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
