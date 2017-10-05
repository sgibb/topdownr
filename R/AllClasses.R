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
#' *Abstract/VIRTUAL* parent class for [TopDownSet-class] and [NCBSet-class] to
#' provide common interface.
#'
#' @details
#' This class just provides a common interface. It is not intended for direct
#' use by the user. Please see [TopDownSet-class] for an example usage of its
#' child class.
#'
#' @slot rowViews [Biostrings::XStringViews-class], information about fragments/bonds
#' (name, type, sequence, mass, charge), see [Biostrings::XStringViews-class] and
#' [FragmentViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class], information about the MS2
#' experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class], intensity/coverage values of the
#' fragments/bonds. The rows correspond to the fragments/bonds and the
#' columns to the condition/run. It just stores values that are different from
#' zero.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`, tolerance in *ppm* that were used for matching the
#' experimental *m/z* values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return This is an *Abstract/VIRTUAL* class to provide a common interface for
#' [TopDownSet-class] and [NCBSet-class]. It is not possible to create an
#' `AbstractTopDownSet` object.
#'
#' @seealso
#' - [TopDownSet-class] and [NCBSet-class] which both implement/use this
#' interface. These manual pages also provide some example code.
#' - [FragmentViews-class] (and [Biostrings::XStringViews-class]) for the row
#' view interface.
#' - [Matrix::dgCMatrix-class] for technical details about the intensity/coverage storage.
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
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
#' See `vignette("analysis", package="topdownr")` for a detailed example how to
#' work with `TopDownSet` objects.
#'
#' @slot rowViews [FragmentViews-class], information about fragments
#' (name, type, sequence, mass, charge), see [FragmentViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class], information about the MS2
#' experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class], intensity values of the fragments. The
#' rows correspond to the fragments and the columns to the condition/run. It
#' just stores values that are different from zero.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`, tolerance in *ppm* that were used for matching the
#' experimental mz values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return An [TopDownSet-class] object.
#'
#' @seealso
#' - [FragmentViews-class] for the row view interface.
#' - [Matrix::dgCMatrix-class] for technical details about the intensity storage.
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
#' @slot rowViews [Biostrings::XStringViews-class], information about bonds
#' (name, start, end, width, sequence), see [Biostrings::XStringViews-class]
#' for details.
#' @slot colData [S4Vectors::DataFrame-class], information about the MS2
#' experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class], coverage values of the bonds. The
#' rows correspond to the bonds and the columns to the condition/run. It
#' just stores values that are different from zero. If a bond is covered by an
#' N-terminal fragment its encoded with `1`, by an C-terminal fragmentl with `2`
#' and by both fragment types by `3` respectively.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`, tolerance in *ppm* that were used for matching the
#' experimental mz values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return An [NCBSet-class] object.
#'
#' @seealso
#' - An `NCBSet` is generated from an [TopDownSet-class] object.
#' - [Biostrings::XStringViews-class] for the row view interface.
#' - [Matrix::dgCMatrix-class] for technical details about the coverage storage.
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
setClass("NCBSet",
         contains="AbstractTopDownSet",
         validity=function(object).validateNCBSet(object)
)
