#' Investigation of Fragmentation Conditions in Top-Down Proteomics
#'
#' The topdownr package allows automatic and systemic investigation of
#' fragment conditions. It creates Thermo Orbitrap Fusion Lumos method files
#' to test hundreds of fragmentation conditions. Additionally it provides
#' functions to analyse and process the generated MS data and determine the
#' best conditions to maximise overall fragment coverage.
#'
#' The usage of the topdownr package is demonstrated in the following vignettes:
#'
#' - Generate .meth files prior data acquisition for the Thermo Orbitrap Fusion
#' Lumos MS devise: `vignette("data-generation", package="topdownr")`.
#' - How to analyse top-down fragmenation data:
#' `vignette("analysis", package="topdownr")`
#'
#' @docType package
#' @name topdownr-package
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de},
#' Pavel Shliaha \email{pavels@bmb.sdu.dk},
#' Ole Nørregaard Jensen \email{jenseno@bmb.sdu.dk}
#' @references \url{https://github.com/sgibb/topdownr/}
#' @keywords package
#'
#' @import methods
#' @import BiocGenerics
#' @import ProtGenerics
#' @import Biostrings
#' @import S4Vectors
#' @importClassesFrom Matrix Matrix dgCMatrix
#' @importFrom Matrix Matrix sparseMatrix sparseVector crossprod tcrossprod
#' drop0 nnzero
#' @importFrom Biobase assayData
#' @importFrom ggplot2 aes aes_string element_blank element_rect element_text
#' facet_grid geom_hline geom_raster geom_segment geom_text geom_vline ggplot
#' ggtitle labs scale_alpha scale_color_manual scale_fill_manual
#' scale_x_discrete scale_y_continuous theme theme_classic
#' @importFrom grDevices pdf dev.off
#' @importFrom mzR openMSfile close header peaks runInfo
#' @importFrom MSnbase calculateFragments defaultNeutralLoss get.amino.acids
#' @importClassesFrom MSnbase MSnSet
#' @importFrom stats ave median setNames quantile
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils modifyList txtProgressBar setTxtProgressBar object.size
#' read.csv .DollarNames packageVersion
"_PACKAGE"
