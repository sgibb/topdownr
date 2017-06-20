#' Internal default function to define default values.
#'
#' @param \ldots named arguments that should be overwritten
#' @param default named list of default arguments
#' @param valid character vector that contains valid tags
#' @return list with default/modified settings
#' @noRd
.defaultSettings <- function(..., default, valid) {
  dots <- list(...)
  m <- match(names(dots), valid)
  if (anyNA(m)) {
    stop("The following setting(s) is/are not valid: ",
         paste0(names(dots)[is.na(m)]))
  }
  modifyList(default, dots)
}

#' Settings for MS1 parameters.
#'
#' @param \ldots named arguments that should be overwritten
#' @return named \code{list} of settings
#' @export
#' @examples
#' library("topdown")
#'
#' # all default settings
#' defaultMs1Settings()
#'
#' # overwrite FirstMass and set AgcTarget
#' defaultMs1Settings(FirstMass=100, AgcTarget=c(1e4, 1e5))
defaultMs1Settings <- function(...) {
  .defaultSettings(...,
                   default=list(FirstMass=400,
                                LastMass=1200,
                                Microscans=10),
                   valid=.validMs1Tags())
}

#' Settings for MS2 parameters.
#'
#' @param \ldots named arguments that should be overwritten
#' @return named \code{list} of settings
#' @export
#' @examples
#' library("topdown")
#'
#' # all default settings
#' defaultMs2Settings()
#'
#' # overwrite AgcTarget
#' defaultMs2Settings(AgcTarget=c(1e4, 1e5))
defaultMs2Settings <- function(...) {
  .defaultSettings(..., default=list(
      ActivationType="ETD",
      AgcTarget=c(1e5, 5e5, 1e6),
      ETDReagentTarget=c(1e6, 5e6, 1e7),
      ETDReactionTime=c(0, 2.5, 5, 10, 15, 30, 50),
      ETDSupplementalActivation=c("ETciD", "EThcD"),
      ETDSupplementalActivationEnergy=seq(0, 35, by=7),
      IsolationWindow=1,
      MaxITTimeInMS=200,
      Microscans=40,
      OrbitrapResolution="R120K"),
    valid=.validMs2Tags())
}
