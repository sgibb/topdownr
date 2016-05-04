defaultProteins <- list(
  myoglobin=matrix(c(738.01, 808.15, 893.16,
                     rep(10, 3)), ncol=2, dimnames=list(c(), c("mass", "z"))),
  h4       =matrix(c(700.5, 737.32, 778.22,
                     rep(10, 3)), ncol=2, dimnames=list(c(), c("mass", "z"))),
  h2a      =matrix(c(609.21, 667.19, 737.26,
                     23, 21, 19), ncol=2, dimnames=list(c(), c("mass", "z"))),
  h2b      =matrix(c(600.51, 657.56, 726.72,
                     23, 21, 19), ncol=2, dimnames=list(c(), c("mass", "z"))),
  h33tail  =matrix(c(446.10, 486.56, 535.11, 594.46, 668.64,
                     12:8), ncol=2, dimnames=list(c(), c("mass", "z")))
)
