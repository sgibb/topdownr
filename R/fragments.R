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
  d <- calculateFragments(sequence,
                          type=type,
                          modifications=modifications,
                          neutralLoss=neutralLoss,
                          verbose=verbose)
  n <- nchar(sequence)
  FragmentViews(sequence, mass=d$mz, type=d$type, z=Rle(d$z), names=d$ion,
                start=ifelse(startsWith(sequence, d$seq), 1L, n - d$pos + 1L),
                width=d$pos)
}
