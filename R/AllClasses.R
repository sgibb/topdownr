#' The FragmentViews class
#'
#' The FragmentViews class is a basic container for storing a set of views
#' (start/end locations) on the same peptides/protein sequence. Additionally
#' it keeps information about mass, type and charge of the fragments.
#'
#' @details FragmentViews extends [Biostrings::XStringViews-class].
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

#' The AbstractTopDownSet class
#'
#' Abstract/VIRTUAL parent class for [TopDownSet-class] and [NCBSet-class] to
#' provide common interface.
#'
#' @noRd
setClass("AbstractTopDownSet",
         contains="VIRTUAL",
         slots=c(rowViews="XStringViews",
                 colData="DataFrame",
                 assay="dgCMatrix",
                 files="character",
                 tolerance="numeric",
                 processing="character"),
         prototype=prototype(rowViews=new("XStringViews"),
                             colData=new("DataFrame"),
                             assay=new("dgCMatrix"),
                             files=character(),
                             tolerance=double(),
                             processing=character()),
         validity=function(object).validateAbstractTopDownSet(object)
)


#' The TopDownSet class
#'
#' The TopDownSet class is a container for a whole top-down proteomics
#' experiment.
#'
#' @details
#' See `vignette("TopDown", package="topdown")` for a detailed example how to
#' work with `TopDownSet` objects.
#'
#' @slot rowViews [FragmentViews-class], information about fragments
#' (name, type, sequence, mass, charge), see [FragmentViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class], information about the MS2
#' experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class], intensity values of the fragments. The
#' rows corresponding to the fragments and the columns to the condition/run. It
#' just stores values that are different from zero.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`, tolerance in *ppm* that were used for matching the
#' experimental mz values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @seealso [FragmentViews-class]
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
setClass("TopDownSet",
         contains="AbstractTopDownSet",
         prototype=prototype(rowViews=new("FragmentViews"))
)

#' The NCBSet class
#'
#' The NCBSet class is a container for a top-down proteomics
#' experiment similar to the [TopDownSet-class] but instead of intensity values
#' it just stores the information if a bond is covered by a
#' N-terminal (encoded as `1`), a C-terminal (encoded as `2`)
#' and/or both fragments (encoded as `3`).
#'
#' @seealso [TopDownSet-class]
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
setClass("NCBSet",
         contains="AbstractTopDownSet",
         validity=function(object).validateNCBSet(object)
)
