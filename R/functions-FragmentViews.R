#' @describeIn FragmentViews-class Constructor
#'
#' In general it is not necessary to call the constructor manually. See
#' [readTopDownFiles()] instead.
#'
#' @param sequence `character`/`AAString`, complete protein/peptide sequence.
#' @param mass `double`, mass of the fragments, same length as
#' `start`/`end`/`width`.
#' @param type `character`, type of the fragments, same length as
#' `start`/`end`/width`.
#' @param z `integer`, charge of the fragments, length one or same length as
#' `start`/`end`/width`.
#' @param start `integer`, start positions of the fragments. At least two of
#' `start`/`end`/width` has to be given.
#' @param end `integer`, end positions of the fragments. At least two of
#' `start`/`end`/width` has to be given.
#' @param width `integer`, width positions of the fragments. At least two of
#' `start`/`end`/width` has to be given.
#' @param names `character`, names of the fragments, same length as
#' `start`/`end`/width`.
#' @param metadata `list`, metadata like modifications.
#' @return An [FragmentViews-class] object.
#' @export
#' @examples
#'
#' # constructor
#' fv <- FragmentViews("ACE", start=1, width=1:3, names=paste0("b", 1:3),
#'                     mass=c(72.04439, 232.07504, 361.11763),
#'                     type="b", z=1)
#' fv
FragmentViews <- function(sequence, mass, type, z=1L,
                          start=NULL, end=NULL, width=NULL, names=NULL,
                          metadata=list()) {
    v <- Views(AAString(sequence), start=start, end=end, width=width,
               names=names)
    d <- DataFrame(mass=mass, type=factor(type), z=Rle(z))
    elementMetadata(v) <- d
    metadata(v) <- metadata
    new("FragmentViews", v[order(d$mass)])
}

#' Validate FragmentViews
#'
#' @param object FragmentViews
#' @return `TRUE` (if valid) else `character` with msg what was incorrect
#' @noRd
.validateFragmentViews <- function(object) {
    msg <- character()

    cols <- c("mass", "type", "z")
    if (!all(cols %in% colnames(elementMetadata(object)))) {
        sel <- !cols %in% colnames(elementMetadata(object))
        msg <- c(msg, paste(paste0("'", cols[sel], "'", collapse=", "),
                            if (sum(sel) == 1L) { "is" } else { "are" },
                            "missing."))
    }

    mass <- elementMetadata(object)[, "mass", drop=TRUE]

    if (!is.double(mass)) {
        msg <- c(msg, "'mass' has to be of type double.")
    }

    if (is.unsorted(mass)) {
        msg <- c(msg, "'mass' has to be sorted.")
    }

    if (length(msg)) {
        msg
    } else {
        TRUE
    }
}
