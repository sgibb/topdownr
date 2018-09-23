#' Find xml children names and tags
#'
#' @param xsd `xml2::xml_document`
#' @param name `character`, element name
#'
#' @return `matrix` of valid tag names and their types
#' @noRd
.xmlChildrenNameType <- function(xsd, name) {
    requireNamespace("xml2")
    nds <- xml2::xml_find_all(
        xsd,
        sprintf(".//xs:element[@name = '%s']//xs:element", name)
    )
    attrs <- xml2::xml_attrs(nds)
    m <- do.call(rbind, lapply(attrs, "[", c("name", "type")))
    repl <- matrix(c("double$", "double",
                     "int.*$", "integer",
                     "boolean$", "logical",
                     "string$", "character"), byrow=TRUE, nrow=4)
    repl[, 1L] <- paste0("^xs:", repl[, 1L])
    for (i in seq_len(nrow(repl))) {
        m[, "type"] <- sub(repl[i, 1L], repl[i, 2L], m[, "type"])
    }
    colnames(m) <- c("name", "class")
    m
}

#' Find xml enumerations
#'
#' @param xsd `xml2::xml_document`
#' @param name `character`, element name
#'
#' @return `character` of valid enumerations
#' @noRd
.xmlEnumeration <- function(xsd, name) {
    requireNamespace("xml2")
    nds <- xml2::xml_find_all(
        xsd,
        sprintf(".//xs:simpleType[@name = '%s']//xs:enumeration", name)
    )
    xml2::xml_attr(nds, "value")
}

#' Find TMSnScan parameters
#'
#' @param xsd `xml2::xml_document`
#'
#' @return `matrix` (`character`)
#' @noRd
.xmlTmsnScanParamters <- function(xsd) {
    param <- .xmlChildrenNameType(xsd, "TMSnScan")
    isEnum <-
        !param[, "class"] %in% c("logical", "integer", "double", "character")
    param[isEnum, "class"] <- vapply(
        param[isEnum, "class"],
        function(tt)paste0(.xmlEnumeration(xsd, tt), collapse=":"),
        NA_character_
    )
    # Exclude "Ms2 parameters"
    param[!startsWith(param[, "name"], "Ms2"),]
}

#' Find valid MS1/MSn settings
#'
#' @param xsd `xml2::xml_document`
#' @return `matrix` (`character`), 3 columns: name, class, type
#' @noRd
.xmlValidMsSettings <- function(xsd) {
    settings <- rbind(
        cbind(.xmlChildrenNameType(xsd, "FullMSScan"), type="MS1"),
        cbind(.xmlTmsnScanParamters(xsd), type="MS2"))

    for (actType in c("ETD", "CID", "HCD", "UVPD")) {
        isType <- grepl(actType, settings[, "name"])
        settings[isType, "type"] <- actType
    }
    settings
}
