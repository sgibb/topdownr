#' Add adducts to the output of MSnbase::calculateFragments
#'
#' @param x data.frame, output of MSnbase::calculateFragments
#' @param adducts data.frame, with 3 columns mass, name, to
#' @return data.frame
#' @noRd
.addAdducts <- function(x, adducts) {
  if (!nrow(adducts)) {
    return(x)
  }

  if (!all(c("mass", "name", "to") %in% colnames(adducts))) {
    stop("The 'adducts' data.frame must have the columns: ",
         "'mass', 'name' and 'to'.")
  }

  r <- do.call(rbind, lapply(seq_len(nrow(adducts)), function(i) {
    a <- x[x$type == adducts$to[i], , drop=FALSE]
    a$mz <- a$mz + adducts$mass[i]
    a$ion <- paste0(adducts$name[i], a$pos)
    a
  }))
  x <- rbind(x, r)
  rownames(x) <- NULL
  x
}

#' Calculate Fragments (via MSnbase::calculateFragments)
#'
#' @param sequence character, peptide sequence
#' @param type character, type of fragments (should be some more than in
#' MSnbase)
#' @param modification modifications (see MSnbase::calculateFragments)
#' @param neutralLoss neutral loss (see MSnbase::calculateFragments)
#' @param adducts data.frame, with 3 columns mass, name, to
#' @param verbose logical, verbose output?
#' @return FragmentViews
#' @noRd
.calculateFragments <- function(sequence, type=c("a", "b", "c", "x", "y", "z"),
                                modifications=c(C=57.02146),
                                neutralLoss=defaultNeutralLoss(),
                                adducts=data.frame(),
                                verbose=interactive()) {
  csequence <- as.character(sequence)
  d <- calculateFragments(csequence,
                          type=type,
                          modifications=modifications,
                          neutralLoss=neutralLoss,
                          verbose=verbose)
  d <- .addAdducts(d, adducts)

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
