#' The FragmentViews class
#'
#' The FragmentViews class is a basic container for storing a set of views
#' (start/end locations) on the same peptides/protein sequence. Additionally
#' it keeps information about mass, type and charge of the fragments.
#'
#' @details FragmentViews extends
#' [Biostrings::XStringViews-class].
#' In short it combines an
#' [IRanges::IRanges-class]
#' object to store start/end
#' location on a sequence, an
#' [Biostrings::AAString] object.
#'
#' @seealso [Biostrings::XStringViews-class]
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' # Constructor
#' fv <- FragmentViews("ACE", start=1, width=1:3, names=paste0("b", 1:3),
#'                     mass=c(72.04439, 232.07504, 361.11763),
#'                     type="b", z=1)
#' fv
#'
#' # Coercion to data.frame
#' as(fv, "data.frame")
## just subclass to ensure that the elementMetadata slot (DataFrame) contains
## a mass column
setClass(
    "FragmentViews",
    contains="XStringViews",
    validity=function(object).validateFragmentViews(object)
)

#' The AbstractTopDownSet class
#'
#' *Abstract/VIRTUAL* parent class for
#' [TopDownSet-class] and [NCBSet-class] to provide common interface.
#'
#' @details
#' This class just provides a common interface. It is not intended for direct
#' use by the user. Please see [TopDownSet-class] for an example usage of its
#' child class.
#'
#' @slot rowViews [Biostrings::XStringViews-class],
#' information about fragments/bonds (name, type, sequence, mass, charge),
#' see [Biostrings::XStringViews-class] and
#' [FragmentViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class],
#' information about the MS2 experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class],
#' intensity/coverage values of the fragments/bonds.
#' The rows correspond to the fragments/bonds and the
#' columns to the condition/run. It just stores values that are
#' different from zero.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`,
#' tolerance in *ppm* that were used for matching the experimental
#' *m/z* values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return This is an *Abstract/VIRTUAL* class
#' to provide a common interface for
#' [TopDownSet-class] and [NCBSet-class].
#' It is not possible to create an `AbstractTopDownSet` object.
#'
#' @seealso
#' - [TopDownSet-class] and [NCBSet-class]
#' which both implement/use this interface.
#' These manual pages also provide some example code.
#' - [FragmentViews-class]
#' (and [Biostrings::XStringViews-class])
#' for the row view interface.
#' - [Matrix::dgCMatrix-class]
#' for technical details about the intensity/coverage storage.
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' ## Because AbstractTopDownSet is a VIRTUAL class we could not create any
#' ## object of it. Here we demonstrate the usage with an TopDownSet that
#' ## implements the AbstractTopDownSet interface. See `?"TopDownSet-class"` for
#' ## more details an further examples.
#'
#' ## Example data
#' data(tds, package="topdownr")
#'
#' tds
#'
#' head(summary(tds))
#'
#' # Accessing slots
#' rowViews(tds)
#' colData(tds)
#' assayData(tds)
#'
#' # Accessing colData
#' tds$Mz
#' tds$FilterString
#'
#' # Subsetting
#'
#' # First 100 fragments
#' tds[1:100]
#'
#' # All c fragments
#' tds["c"]
#'
#' # Just c 152
#' tds["c152"]
#'
#' # Condition 1 to 10
#' tds[, 1:10]
#
setClass(
    "AbstractTopDownSet",
    contains="VIRTUAL",
    slots=c(
        rowViews="XStringViews",
        colData="DataFrame",
        assay="dgCMatrix",
        files="character",
        tolerance="numeric",
        processing="character"
    ),
    prototype=prototype(
        rowViews=new("XStringViews"),
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
#' See `vignette("analysis", package="topdownr")`
#' for a detailed example how to work with
#' `TopDownSet` objects.
#'
#' @slot rowViews [FragmentViews-class],
#' information about fragments (name, type, sequence, mass, charge), see
#' [FragmentViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class],
#' information about the MS2 experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class],
#' intensity values of the fragments. The
#' rows correspond to the fragments and the columns to the condition/run. It
#' just stores values that are different from zero.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`,
#' tolerance in *ppm* that were used for matching the
#' experimental mz values to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return An [TopDownSet-class] object.
#'
#' @seealso
#' - [FragmentViews-class]
#' for the row view interface.
#' - [Matrix::dgCMatrix-class]
#' for technical details about the intensity storage.
#' - `?vignette("analysis", package="topdownr")`
#' for a full documented example of an analysis using `topdownr`.
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' ## Example data
#' data(tds, package="topdownr")
#'
#' tds
#'
#' head(summary(tds))
#'
#' # Accessing slots
#' rowViews(tds)
#' colData(tds)
#' assayData(tds)
#'
#' # Accessing colData
#' tds$Mz
#' tds$FilterString
#'
#' # Subsetting
#'
#' # First 100 fragments
#' tds[1:100]
#'
#' # All c fragments
#' tds["c"]
#'
#' # Just c 152
#' tds["c152"]
#'
#' # Condition 1 to 10
#' tds[, 1:10]
#'
#' # Filtering
#' # Filter all intensities that don't have at least 10 % of the highest
#' # intensity per fragment.
#' tds <- filterIntensity(tds, threshold=0.1)
#'
#' # Filter all conditions with a CV above 30 % (across technical replicates)
#' tds <- filterCv(tds, threshold=30)
#'
#' # Filter all conditions with a large deviation in injection time
#' tds <- filterInjectionTime(tds, maxDeviation=log2(3), keepTopN=2)
#'
#' # Filter all conditions where fragments don't replicate
#' tds <- filterNonReplicatedFragments(tds)
#'
#' # Normalise by TIC
#' tds <- normalize(tds)
#'
#' # Aggregate technical replicates
#' tds <- aggregate(tds)
#'
#' head(summary(tds))
#'
#' # Coercion
#' as(tds, "NCBSet")
#'
#' if (require("MSnbase")) {
#'     as(tds, "MSnSet")
#' }
setClass(
    "TopDownSet",
    contains="AbstractTopDownSet",
    prototype=prototype(rowViews=new("FragmentViews"))
)

#' The NCBSet class
#'
#' The NCBSet class is a container for a top-down proteomics experiment
#' similar to the [TopDownSet-class]
#' but instead of intensity values it just stores the information if a
#' bond is covered by a N-terminal (encoded as `1`),
#' a C-terminal (encoded as `2`)
#' and/or both fragments (encoded as `3`).
#'
#' @slot rowViews [Biostrings::XStringViews-class],
#' information about bonds (name, start, end, width, sequence),
#' see [Biostrings::XStringViews-class] for details.
#' @slot colData [S4Vectors::DataFrame-class],
#' information about the MS2 experiments and the fragmentation conditions.
#' @slot assay [Matrix::dgCMatrix-class],
#' coverage values of the bonds. The rows correspond to the bonds and the
#' columns to the condition/run. It just stores values that are different
#' from zero. If a bond is covered by an N-terminal fragment its encoded
#' with `1`, by an C-terminal fragmentl with `2` and
#' by both fragment types by `3` respectively.
#' @slot files `character`, files that were imported.
#' @slot tolerance `double`,
#' tolerance in *ppm* that were used for matching the experimental mz values
#' to the theoretical fragments.
#' @slot processing `character`, log messages.
#'
#' @return An [NCBSet-class] object.
#'
#' @seealso
#' - An `NCBSet` is generated from an [TopDownSet-class] object.
#' - [Biostrings::XStringViews-class]
#' for the row view interface.
#' - [Matrix::dgCMatrix-class]
#' for technical details about the coverage storage.
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' ## Example data
#' data(tds, package="topdownr")
#'
#' ## Aggregate technical replicates
#' tds <- aggregate(tds)
#'
#' ## Coercion into an NCBSet object
#' ncb <- as(tds, "NCBSet")
#'
#' ncb
#'
#' head(summary(ncb))
#'
#' # Accessing slots
#' rowViews(ncb)
#' colData(ncb)
#' assayData(ncb)
#'
#' # Accessing colData
#' ncb$Mz
#'
#' # Subsetting
#'
#' # First 100 bonds
#' ncb[1:100]
#'
#' # Just bond 152
#' ncb["bond152"]
#'
#' # Condition 1 to 10
#' ncb[, 1:10]
#'
#' # Plot fragmentation map
#' fragmentationMap(ncb)
setClass(
    "NCBSet",
    contains="AbstractTopDownSet",
    validity=function(object).validateNCBSet(object)
)
