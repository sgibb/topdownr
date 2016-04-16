defaultMs1Settings <- function(...) {
  default <- list(FirstMass=550, LastMass=1200,
                  Microscans=10)
  validTags <- c("FirstMass", "LastMass", "Microscans",
                 "MaxITTimeInMS", "AgcTarget")
  dots <- list(...)
  dots <- dots[names(dots) %in% validTags]
  modifyList(default, dots)
}
