#' helper function to create default proteins entries
#'
#' @param x double, mass values for the given protein
#' @return named matrix length(x)x2
#' @noRd
.mzMatrix <- function(x) {
  matrix(c(x, rep(10L,
                  # regardless of the true charge Xcalibur expects 10 here
                  length(x))),
         ncol=2L, dimnames=list(c(), c("mass", "z")))
}

#' list of default proteins used to created default method files easier
#' @noRd
defaultProteins <- list(
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


