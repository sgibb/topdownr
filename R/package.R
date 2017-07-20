#' Helper functions for topdown proteomics.
#'
#' Some xml helper functions for topdown proteomics.
#'
#' @docType package
#' @name topdown-package
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @references \url{https://github.com/sgibb/topdown/}
#' @keywords package
#'
#' @import methods
#' @import S4Vectors
#' @import Biostrings
#' @importClassesFrom Matrix Matrix dgCMatrix
#' @importFrom Matrix Matrix sparseMatrix sparseVector nnzero drop0
#' @importFrom Biobase assayData
#'
#' @importFrom mzR openMSfile close header peaks
#' @importFrom MSnbase calculateFragments defaultNeutralLoss
#' @importFrom stats setNames
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils modifyList txtProgressBar setTxtProgressBar object.size read.csv .DollarNames
NULL
