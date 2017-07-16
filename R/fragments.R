#' Calculate Fragments (via MSnbase::calculateFragments)
#'
#' @param sequence character, peptide sequence
#' @param type character, type of fragments (should be some more than in
#' MSnbase)
#' @param modification modifications (see MSnbase::calculateFragments)
#' @param neutralLoss neutral loss (see MSnbase::calculateFragments)
#' @param verbose logical, verbose output?
#' @return FragmentViews
#' @noRd
.calculateFragments <- function(sequence, type=c("a", "b", "c", "x", "y", "z"),
                                modifications=c(C=57.02146),
                                neutralLoss=defaultNeutralLoss(),
                                verbose=interactive()) {
  csequence <- as.character(sequence)
  d <- calculateFragments(csequence,
                          type=type,
                          modifications=modifications,
                          neutralLoss=neutralLoss,
                          verbose=verbose)
  n <- nchar(sequence)
  FragmentViews(sequence, mass=d$mz, type=d$type, z=Rle(d$z), names=d$ion,
                start=ifelse(startsWith(csequence, d$seq), 1L, n - d$pos + 1L),
                width=d$pos)
}

#' Match fragments and measured mz values.
#'
#' @param mz double, measured mz.
#' @param fmass double, fragment mass
#' @param tolerance double, allowed tolerance
#' @return integer
#' @noRd
.matchFragments <- function(mz, fmass, tolerance=2.5e-5) {
  if (!length(mz)) {
    integer()
  }
  m <- MSnbase:::relaxedMatch(mz, fmass, nomatch=NA, tolerance=tolerance,
                              relative=TRUE)
  if (anyDuplicated(m)) {
    o <- order(abs(mz - fmass[m]))
    sortedMatches <- m[o]
    sortedMatches[which(duplicated(sortedMatches))] <- NA
    m[o] <- sortedMatches
  }
  as.integer(m)
}
