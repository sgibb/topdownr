#' Find xml children names and tags
#'
#' @return `data.frame`, of valid tag names and their types
#' @noRd
.xmlChildrenNameType <- function(xsd, name) {
    nds <- xml2::xml_find_all(
        xsd,
        sprintf(".//xs:element[@name = '%s']//xs:element", name)
    )
    attrs <- xml2::xml_attrs(nds)
    m <- do.call(rbind, lapply(attrs, "[", c("name", "type")))
    m[, "type"] <- gsub("xs:double", "double", m[, "type"], fixed=TRUE)
    m[, "type"] <- gsub("xs:int", "integer", m[, "type"], fixed=TRUE)
    m[, "type"] <- gsub("xs:boolean", "logical", m[, "type"], fixed=TRUE)
    m
}

#' Find xml enumerations
#'
#' @return `character` of valid enumerations
#' @noRd
.xmlEnumeration <- function(xsd, name) {
    nds <- xml2::xml_find_all(
        xsd,
        sprintf(".//xs:simpleType[@name = '%s']//xs:enumeration", name)
    )
    xml2::xml_attr(nds, "value")
}
