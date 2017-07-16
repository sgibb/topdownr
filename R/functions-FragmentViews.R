#' @export
FragmentViews <- function(sequence, mass, type, z=1L,
                          start=NULL, end=NULL, width=NULL, names=NULL) {
  v <- Views(AAString(sequence), start=start, end=end, width=width, names=names)
  d <- DataFrame(mass=mass, type=factor(type), z=Rle(z))
  elementMetadata(v) <- d
  new("FragmentViews", v[order(d$mass)])
}

#' Validate FragmentViews
#'
#' @param object FragmentViews
#' @return TRUE (if valid) else character with msg what was incorrect
#' @noRd
.validateFragmentViews <- function(object) {
  msg <- character()

  cols <- c("mass", "type", "z")
  if (!all(cols %in% colnames(elementMetadata(object)))) {
    sel <- !cols %in% colnames(elementMetadata(object))
    msg <- c(msg, paste(paste0("'", cols[sel], "'", collapse=", "),
                        if (sum(sel) == 1L) { "is" } else { "are" },
                        "missing."))
  }

  mass <- elementMetadata(object)[, "mass", drop=TRUE]

  if (!is.double(mass)) {
    msg <- c(msg, "'mass' has to be of type double.")
  }

  if (is.unsorted(mass)) {
    msg <- c(msg, "'mass' has to be sorted.")
  }

  if (length(msg)) {
    msg
  } else {
    TRUE
  }
}
