.defaultSettings <- function(..., default, validFun) {
  validFun <- match.fun(validFun)
  dots <- list(...)
  dots <- dots[names(dots) %in% validFun()]
  modifyList(default, dots)
}

defaultMs1Settings <- function(...) {
  .defaultSettings(..., default=list(FirstMass=550, LastMass=1200,
                                     Microscans=10),
                   validFun=.validMs1Tags)
}

defaultMs2Settings <- function(...) {
  .defaultSettings(..., default=list(IsolationWindow=5, MaxITTimeInMS=150,
                                     OrbitrapResolution="R60K",
                                     Microscans=20),
                   validFun=.validMs2Tags)
}

.validMs1Tags <- function() {
  c("FirstMass", "LastMass", "Microscans",
    "MaxITTimeInMS", "AgcTarget")
}

.validMs2Tags <- function() {
  c("ActivationType", "IsolationWindow", "EnableMultiplexIons", "EnableMSXIds",
    "MaxNoOfMultiplexIons", "OrbitrapResolution", "AgcTarget", "MinAgcTarget",
    "MaxITTimeInMS", "Microscans", "ETDReactionTime", "ETDReagentTarget",
    "MaximumETDReagentInjectionTime", "UseInternalCalibratedETD",
    "ETDSupplementalActivationEnergy", "ETDSupplementalActivation")
}

.commonSettings <- function(x) {
  x[lengths(x) == 1L]
}

.nonCommonSettings <- function(x) {
  x[lengths(x) > 1L]
}
