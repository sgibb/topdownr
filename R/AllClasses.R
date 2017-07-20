#' The FragmentViews class
#'
#' The FragmentViews class is a basic container for storing a set of views
#' (start/end locations) on the same peptides/protein sequence. Additionally
#' it keeps information about mass, type and charge of the fragments.
#'
#' @details FragmentViews is inherits [Biostrings::XStringViews-class].
#' In short it combines an [IRanges::IRanges-class] object to store start/end
#' location on a sequence, an [Biostrings::AAString] object.
#'
#' @seealso [Biostrings::XStringViews-class]
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
## just subclass to ensure that the elementMetadata slot (DataFrame) contains
## a mass column
setClass("FragmentViews",
         contains="XStringViews",
         validity=function(object).validateFragmentViews(object))

setClass("TopDownSet",
  slots=c(rowViews="FragmentViews",
          colData="DataFrame",
          assay="Matrix",
          files="character",
          tolerance="numeric",
          processing="character"),
  prototype=prototype(rowViews=new("FragmentViews"),
                      colData=new("DataFrame"),
                      assay=new("dgCMatrix"),
                      files=character(),
                      tolerance=double(),
                      processing=character()),
  validity=function(object).validateTopDownSet(object)
)
