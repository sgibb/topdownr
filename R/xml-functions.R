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
        sprintf(".//xs:element[@name = '%s']//xs:element[@name|@ref]", name)
    )
    attrs <- lapply(xml2::xml_attrs(nds), function(a) {
        name <- if (is.na(a["name"])) a["ref"] else a["name"]
        type <- a["type"]
        if (is.na(type) && name == "MassList") {
            type <- "MassListRecord"
        } else if (is.na(type) && name == "ScanDescription") {
            type <- "character"
        }
        c(name, type)
    })
    m <- do.call(rbind, attrs)
    colnames(m) <- c("name", "class")
    repl <- matrix(c("double$", "double",
                     "int.*$", "integer",
                     "boolean$", "logical",
                     "string$", "character"), byrow=TRUE, nrow=4L)
    repl[, 1L] <- paste0("^xs:", repl[, 1L])
    for (i in seq_len(nrow(repl))) {
        m[, "class"] <- sub(repl[i, 1L], repl[i, 2L], m[, "class"])
    }
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
.xmlTmsnScanParameters <- function(xsd) {
    param <- .xmlChildrenNameType(xsd, "TMSnScan")
    isEnum <- !param[, "class"] %in%
        c("logical", "integer", "double", "character", "MassListRecord")
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
        cbind(.xmlTmsnScanParameters(xsd), type="TMS2"))

    for (actType in c("ETD", "CID", "HCD", "UVPD")) {
        isType <- grepl(actType, settings[, "name"])
        settings[isType, "type"] <- actType
    }
    settings
}
