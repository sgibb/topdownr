#' Calculate Fragments (via MSnbase::calculateFragments)
#'
#' @param sequence character, peptide sequence
#' @param type character, type of fragments (should be some more than in
#' MSnbase)
#' @param modification modifications (see MSnbase::calculateFragments)
#' @param neutralLoss neutral loss (see MSnbase::calculateFragments)
#' @param verbose logical, verbose output?
#' @return data.table
#' @noRd
.calculateFragments <- function(sequence, type=c("a", "b", "c", "x", "y", "z"),
                                modifications=c(C=57.02146),
                                neutralLoss=defaultNeutralLoss(),
                                verbose=interactive()) {
  d <- setDT(calculateFragments(sequence,
                                type=type,
                                modifications=modifications,
                                neutralLoss=neutralLoss,
                                verbose=verbose))
  setorder(d, mz)
  d[, FragmentId := as.double(.I)]
  setkey(d, FragmentId)
  d
}
