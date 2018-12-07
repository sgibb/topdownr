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
        stop(
            "The following setting(s) is/are not valid: ",
            paste0(names(dots)[is.na(m)])
        )
    }
    modifyList(default, dots)
}

#' Settings for MS1/2 parameters.
#'
#' These functions are deprecated. Use [createExperimentsFragmentOptimisation()]
#' instead.
#'
#' @param \ldots named arguments that should be overwritten.
#' @return A named `list` of settings.
#' @rdname default-functions
#' @export
defaultMs1Settings <- function(...) {
    .Deprecated("createExperimentsFragmentOptimisation")
    .defaultSettings(...,
        default=list(
            FirstMass=400,
            LastMass=1200,
            Microscans=10),
        valid=.validMs1Tags()
    )
}

#' @rdname default-functions
#'
#' @export
defaultMs2Settings <- function(...) {
    .Deprecated("createExperimentsFragmentOptimisation")
    .defaultSettings(...,
        default=list(
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
        valid=.validMs2Tags()
    )
}
