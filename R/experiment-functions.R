#' Collapse MassList
#'
#' @param x `matrix`, 2 columns: mz, z
#' @return `character`
#' @noRd
.collapseMassList <- function(x) {
    stopifnot(is.matrix(x))
    stopifnot(ncol(x) == 2L)
    paste0(apply(x, 1L, paste0, collapse="/"), collapse=" ")
}

#' Expand MS1 Conditions
#'
#' TODO:
#'
#' @param type `character`, *ActivationType* for MS2 one of CID, HCD, ETD, or
#' UVPD.
#' @param ... further named arguments, used to create the combination of
#' conditions.
#' @return `data.frame`
#' @export
expandMs1Conditions <- function(...) {
    settings <- .flatten(list(...))
    .validateMsSettings(type="MS1", settings)
    expand.grid(settings, stringsAsFactors=FALSE)
}

#' Expand MS2 Conditions
#'
#' TODO:
#'
#' @param massList `matrix`, 2 columns (mass, z).
#' @param ActivationType `character`, *ActivationType* for MS2 one of CID, HCD, ETD, or
#' UVPD.
#' @param ... further named arguments, used to create the combination of
#' conditions.
#' @return `data.frame`
#' @export
#' @examples
#' expandMs2Conditions(
#'      massList=cbind(mz=c(560.6, 700.5, 933.7), z=rep(1, 3)),
#'      ActivationType="CID",
#'      OrbitrapResolution="R120K",
#'      IsolationWindow=1,
#'      MaxITTimeInMS=200,
#'      Microscans=as.integer(40),
#'      AgcTarget=c(1e5, 5e5, 1e6),
#'      CIDCollisionEnergy=c(NA, seq(7, 35, 7))
#' )
expandMs2Conditions <- function(massList,
                                ActivationType=c("CID", "HCD", "ETD", "UVPD"),
                                ...) {
    ActivationType <- match.arg(ActivationType)
    settings <- .flatten(list(...))
    .validateMsSettings(type=ActivationType, settings)
    expand.grid(
        c(
            MassList=.collapseMassList(massList),
            ActivationType=ActivationType,
            settings
        ),
        stringsAsFactors=FALSE
    )
}

#' Validate a single MS setting against internal .validMsSettings (derivated
#' from XSD)
#'
#' @param name `character`, element name
#' @param value any type, value of element
#' @param type `character`, type of setting
#' @return `TRUE` if valid, else message
#' @noRd
.validateMsSetting <- function(name, value, type) {
    isValidEntry <- .validMsSettings[, "name"] == name &
        (.validMsSettings[, "type"] == type |
         .validMsSettings[, "type"] == "MS2" &
         type != "MS1")

    if (!any(isValidEntry)) {
        return(
            paste0(
                name, " is not a valid element of type '", type, "'",
                " and/or not defined in MethodModification.xsd."
            )
        )
    }

    entry <- .validMsSettings[isValidEntry, "class"]
    cl <- typeof(value)

    if (isTRUE(cl == "character")) {
        isValidValue <-
            .vapply1l(value, function(x)grepl(paste0("(^|:)", x, "(:|$)"), entry))

        if (any(!isValidValue)) {
            return(
                paste0(
                    name, " could not be '", paste0(value[!isValidValue]),
                    "'. Should be one of '", gsub(":", ", ", entry), "'."
                )
            )
        }
    } else if (cl != entry) {
        return(
            paste0(name, " should be of class '", entry, "' but is '", cl, "'.")
        )
    }
    TRUE
}

#' Validate MS settings
#'
#' @param type `character`, MS1/MS2/ActivationType
#' @param settings `list`, named arguments used for validation
#' @return `TRUE` if valid, else stops with an error
#' @noRd
.validateMsSettings <- function(type=c("MS1", "MS2", "ETD", "CID", "HCD", "UVPD"),
                                settings) {
    type <- match.arg(type)
    validation <- mapply(
        .validateMsSetting,
        name=names(settings),
        value=settings,
        type=type,
        SIMPLIFY=FALSE
    )
    isValid <- .vapply1l(validation, isTRUE)

    if (any(!isValid)) {
        stop(paste0(unlist(validation[!isValid]), collapse="\n  "))
    }
    TRUE
}


