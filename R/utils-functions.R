#' cat0, cat with sep="", similar to paste0
#' @param \ldots arguments passed to cat
#' @noRd
cat0 <- function(...) {
  cat(..., sep="", append=TRUE)
}
