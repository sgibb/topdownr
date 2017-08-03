#' Read top-down files.
#'
#' It creates an [TopDownSet-class] object and is its only constructor.
#'
#' @details
#'
#' `readTopDownFiles` reads and processes all top-down files, namely:
#' * `.fasta` (protein sequence)
#' * `.mzML` (spectra)
#' * `.experiments.csv` (method/fragmentation conditions)
#' * `.txt` (scan header information)
#'
#' `adducts`: *Thermo's Xtract* allows some mistakes in deisotoping, mostly it
#' allows `+/- C13-C12` and `+/- H+`. The `adducts` argument takes a
#' `data.frame` with the `mass` to add, the `name` that should assign to these
#' new fragments and an information `to` whom the modification should be
#' applied, e.g. for `H+` on `z`, `data.frame(mass=1.008, name="zpH", to="z")`.
#'
#' *Please note:* The `adducts` are added to the output of
#' [MSnbase::calculateFragments()]. That has some limitations, e.g.
#' neutral loss calculation could not be done in [topdown-package]. If neutral
#' loss should be applied on adducts you have to create additional rows, e.g.:
#' `data.frame(mass=c(1.008, 1.008), name=c("cpH", "cpH_"), to=c("c", "c_"))`.
#'
#' @param path `character`, path to directory that contains the top-down files.
#' @param pattern `character`, a filename pattern, the default `.*` means all
#' files.
#' @param type `character`, type of fragments, currently *a-c* and *x-z* are
#' supported, see [MSnbase::calculateFragments()] for details.
#' @param modifications `character`, unimod names of modifications that should
#' be applied. Currenlty just *Acetyl* (Unimod:1 but just protein N-term),
#' *Carbamidomethyl* (Unimod:4), *Met-loss* (Unimod:765), and
#' *Met-loss+Acetyl* (Unimod:766) are supported. Please note that `c("Acetyl",
#' "Met-loss")` is not identical to `c("Met-loss+Acetyl")`. The former removes
#' M (if followed by A, C, G, P, S, T, or V) and does Acetylation on the
#' protein N-term while the latter just acytelate the N-term if M was followed
#' by A, G, S, T (see
#' http://www.unimod.org/modifications_view.php?editid1=1,
#' http://www.unimod.org/modifications_view.php?editid1=765, and
#' http://www.unimod.org/modifications_view.php?editid1=766 for details).
#' Use `NULL` to disable all modifications.
#' @param adducts `data.frame`, with 3 columns, namely: mass, name, to, see
#' details section.
#' @param neutralLoss `list`, neutral loss that should be applied, see
#' [MSnbase::calculateFragments()] and [MSnbase::defaultNeutralLoss()] for
#' details.
#' @param tolerance `double`, tolerance in *ppm* that is used to match the
#' theoretical fragments with the observed ones.
#' @param dropNonInformativeColumns logical, should columns with just one
#' identical value across all runs be removed?
#' @param sampleColumns `character`, column names of the [colData()] used to
#' define a sample (technical replicate). This is used to add the `Sample`
#' column (used for easier aggregation, etc.).
#' @param verbose `logical`, verbose output?
#' @return A `TopDownSet` object.
#' @export
#' @seealso [MSnbase::calculateFragments()], [MSnbase::defaultNeutralLoss()]
#' @examples
#' library("topdown")
#'
#' if (require("topdowndata")) {
#'   # add H+ to z and no neutral loss of water
#'   tds <- readTopDownFiles(topdowndata::topDownDataPath("myoglobin"),
#'                           adducts=data.frame(mass=1.008, name="zpH", to="z"),
#'                           neutralLoss=MSnbase::defaultNeutralLoss(
#'                            disableWaterLoss=c("Cterm", "D", "E", "S", "T")),
#'                           tolerance=25e-6)
#' }
readTopDownFiles <- function(path, pattern=".*",
                             type=c("a", "b", "c", "x", "y", "z"),
                             modifications=c("Carbamidomethyl",
                                             "Met-loss+Acetyl"),
                             adducts=data.frame(),
                             neutralLoss=MSnbase::defaultNeutralLoss(),
                             tolerance=10e-6,
                             dropNonInformativeColumns=TRUE,
                             sampleColumns=c("Mz", "AgcTarget",
                                             "EtdReagentTarget",
                                             "EtdActivation",
                                             "CidActivation",
                                             "HcdActivation",
                                             "SupplementalActivation",
                                             "SupplementalActivationCe"),
                             verbose=interactive()) {

   files <- .listTopDownFiles(path, pattern=pattern)

   sequence <- .readFasta(files$fasta, verbose=verbose)

   fragmentViews <- .calculateFragments(sequence=sequence,
                                        type=type,
                                        modifications=modifications,
                                        neutralLoss=neutralLoss,
                                        adducts=adducts,
                                        verbose=verbose)

   scanConditions <- do.call(rbind, lapply(files$txt, .readScanHeadsTable,
                                           verbose=verbose))

   headerInformation <- do.call(rbind, lapply(files$csv, .readExperimentCsv,
                                              verbose=verbose))

   mzml <- mapply(.readMzMl,
                  file=files$mzML,
                  scans=split(scanConditions$Scan, scanConditions$File),
                  MoreArgs=list(fmass=elementMetadata(fragmentViews)$mass,
                                tolerance=tolerance, verbose=verbose),
                  SIMPLIFY=FALSE, USE.NAMES=FALSE)

   mzmlHeader <- do.call(rbind, lapply(mzml, "[[", "hd"))

   scanHeadsman <- .mergeScanConditionAndHeaderInformation(scanConditions,
                                                           headerInformation)

   header <- .mergeSpectraAndHeaderInformation(mzmlHeader, scanHeadsman)
   header$Sample <- .groupId(header, cols=sampleColumns)
   header$MedianIonInjectionTimeMs <- .medianIonInjectionTime(header)

   if (dropNonInformativeColumns) {
       header <- .dropNonInformativeColumns(header)
   }

   assay <- do.call(cbind, lapply(mzml, "[[", "m"))
   dimnames(assay) <- list(names(fragmentViews),
                           rownames(header))

   new("TopDownSet",
       rowViews=fragmentViews,
       colData=.colsToRle(as(header, "DataFrame")),
       assay=assay,
       files=basename(unlist(unname(files))),
       tolerance=tolerance,
       processing=.logmsg("Data loaded (tolerance: ",
                          round(tolerance/1e-6, 1L), " ppm)."))
}

#' Test for TopDownSet class
#'
#' @param object object to test
#' @return `TRUE` if object is a TopDownSet otherwise fails with an error
#' @noRd
.isTopDownSet <- function(object) {
    if (!isTRUE(is(object, "TopDownSet"))) {
        stop("'object' has to be an 'TopDownSet' object.")
    }
    TRUE
}

#' Get fragment mass
#'
#' @param object `TopDownSet`
#' @return `double`
#' @noRd
fragmentMass <- function(object) {
    .isTopDownSet(object)
    elementMetadata(object@rowViews)$mass
}

#' Get fragment names
#'
#' @param object `TopDownSet`
#' @return `character`
#' @noRd
fragmentNames <- function(object) {
    .isTopDownSet(object)
    names(object@rowViews)
}

#' Get fragment types
#'
#' @param object `TopDownSet`
#' @return `character`
#' @noRd
fragmentType <- function(object) {
    .isTopDownSet(object)
    elementMetadata(object@rowViews)$type
}

#' Create NCB Map (N-/C-terminal, or both)
#'
#' @param object `TopDownSet`
#' @return `Matrix`, Nterm == 1, Cterm == 2, both == 3
#' @noRd
.ncbMap <- function(object, nterm=c("a", "b", "c"), cterm=c("x", "y", "z")) {
    .isTopDownSet(object)

    w <- width(object@rowViews)
    mn <- mc <- object@assay
    selN <- fragmentType(object) %in% nterm
    selC <- fragmentType(object) %in% cterm
    mn[!selN, ] <- 0L
    mc[!selC, ] <- 0L

    mn <- as(.colSumsGroup(mn, w) > 0L, "dgCMatrix")
    mc <- as(.colSumsGroup(mc, max(w) + 1L - w) > 0L, "dgCMatrix")
    mc@x[] <- 2
    mn + mc
}

#' Add log message.
#'
#' @param object `TopDownSet`
#' @return `TopDownSet`
#' @noRd
.tdsLogMsg <- function(object, ...) {
   .isTopDownSet(object)
   object@processing <- c(object@processing, .logmsg(...))
   object
}

#' Validate `TopDownSet`
#'
#' @param object `TopDownSet`
#' @return `TRUE` (if valid) else character with msg what was incorrect
#' @noRd
.validateTopDownSet <- function(object) {
    msg <- character()

    if (nrow(object@assay) != length(object@rowViews)) {
        msg <- c(msg,
                 "Mismatch between fragment data in 'rowViews' and 'assay'.")
    }

    if (any(rownames(object@assay) != names(object@rowViews))) {
        msg <- c(msg,
                 "Mismatch between fragment names in 'rowViews' and 'assay'.")
    }

    if (ncol(object@assay) != nrow(object@colData)) {
        msg <- c(msg,
                 "Mismatch between condition data in 'colData' and 'assay'.")
    }

    if (any(colnames(object@assay) != rownames(object@colData))) {
        msg <- c(msg,
                 "Mismatch between condition names in 'colData' and 'assay'.")
    }

    if (length(msg)) {
        msg
    } else {
        TRUE
    }
}
