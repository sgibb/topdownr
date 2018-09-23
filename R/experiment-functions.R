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
