#' Internal default function to define default values.
#'
#' @param \ldots named arguments that should be overwritten.
#' @param default `list`, named, default arguments.
#' @param valid `character`, containing valid tags.
#' @return `list` with default/modified settings
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

#' Internal function to create default proteins entries.
#'
#' @param x `double`, mass values for the given protein.
#' @return named `matrix` (`length(x)` times `2`).
#' @noRd
.mzMatrix <- function(x) {
  matrix(c(x, rep(10L,
                  # regardless of the true charge Xcalibur expects 10 here
                  length(x))),
         ncol=2L, dimnames=list(c(), c("mass", "z")))
}

#' Settings for MS1/2 parameters.
#'
#' This functions create the default settings for [writeMethodXmls()].
#'
#' @param \ldots named arguments that should be overwritten.
#' @return A named `list` of settings.
#' @rdname default-functions
#' @export
#' @examples
#' library("topdown")
#'
#' # all default MS1 settings
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

#' @rdname default-functions
#'
#' @export
#' @examples
#'
#' # all default MS2 settings
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

#' @rdname default-functions
#'
#' @param protein `character`, protein name.
#'
#' @export
#' @examples
#'
#' # default mass/z values for H2A
#' defaultProteins("h2a")
defaultProteins <- function(protein=c("GST", "myoglobin",
                                      "h2a", "h2b", "h3_1", "h3_3", "h4",
                                      "c345c_c3", "h33tail")) {
  switch(match.arg(protein),
    "GST"      =.mzMatrix(c(799.74, 906.20, 1045.46)),  # z: 34, 30, 26
    "myoglobin"=.mzMatrix(c(738.01, 808.15, 893.16)),
    "h2a"      =.mzMatrix(c(609.21, 700.45, 823.95)),   # z: 23, 20, 17
    "h2b"      =.mzMatrix(c(600.51, 690.39, 812.10)),   # z: 23, 20, 17
    "h3_1"     =.mzMatrix(c(611.87, 664.98, 728.17)),   # z: 25, 23, 21
    "h3_3"     =.mzMatrix(c(608.87, 691.76, 800.77)),   # z: 25, 22, 19
    "h4"       =.mzMatrix(c(625.19, 750.03, 937.29)),   # z: 18, 15, 12
    "c345c_c3" =.mzMatrix(c(950.51, 1062.22, 1203.78)), # z: 19, 17, 15
    "h33tail"  =.mzMatrix(c(446.10, 486.56, 535.11, 594.46, 668.64)) # z: 12:8
  )
}
