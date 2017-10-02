#' Investigation of Fragmentation Conditions in Top-Down Proteomics
#'
#' The topdown package allows automatic and systemic investigation of
#' fragment conditions. It creates Thermo Orbitrap Fusion Lumos method files
#' to test hundreds of fragmentation conditions. Additionally it provides
#' functions to analyse and process the generated MS data and determine the
#' best conditions to maximise overall fragment coverage.
#'
#' The usage of the topdown package is demonstrated in the following vignettes:
#'
#' - Generate .meth files prior data acquisition for the Thermo Orbitrap Fusion
#' Lumos MS devise: `vignette("data-generation", package="topdown")`.
#' - How to analyse top-down fragmenation data: `vignette("analysis", package="topdown")`
#'
#' @docType package
#' @name topdown-package
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de},
#' Pavel Shliaha \email{pavels@bmb.sdu.dk},
#' Ole Nørregaard Jensen \email{jenseno@bmb.sdu.dk}
#' @references \url{https://github.com/sgibb/topdown/}
#' @keywords package
#'
#' @import methods
#' @import BiocGenerics
#' @import Biostrings
#' @import S4Vectors
#' @importClassesFrom Matrix Matrix dgCMatrix
#' @importFrom Matrix Matrix sparseMatrix sparseVector crossprod tcrossprod
#' drop0 nnzero
#' @importFrom Biobase assayData
#' @importFrom ggplot2 ggplot geom_raster aes_string facet_grid
#' scale_fill_manual scale_x_discrete scale_y_continuous geom_vline geom_hline
#' ggtitle theme element_text element_blank element_rect
#'
#' @importFrom mzR openMSfile close header peaks
#' @importFrom MSnbase calculateFragments defaultNeutralLoss get.amino.acids
#' @importClassesFrom MSnbase MSnSet
#' @importFrom stats ave median setNames quantile
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils modifyList txtProgressBar setTxtProgressBar object.size
#' read.csv .DollarNames
"_PACKAGE"
