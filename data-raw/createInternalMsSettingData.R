devtools::load_all()

baseurl <-
    "https://raw.githubusercontent.com/thermofisherlsms/meth-modifications"

## Currently we just support the 3.x serie of the "Calcium" family
version <- list(
    Calcium=c("3.1", "3.2", "3.3")
)

## 3.2/3.3 are (internal) beta (3.1 is released (public) beta?)
release <- list(
    Calcium=paste(c("release/beta", "beta", "beta"), version$Calcium, sep="-")
)

.validMsSettingsXsd <- vector(mode="list", length=length(version))
names(.validMsSettingsXsd) <- names(version)

for (i in seq(along=.validMsSettingsXsd)) {
    .validMsSettingsXsd[[i]] <- vector(mode="list", length=length(version[[i]]))
    names(.validMsSettingsXsd[[i]]) <- version[[i]]

    for (j in seq(along=version[[i]])) {
        url <- paste(
            baseurl,
            release[[i]][j],
            "xsds",
            names(version)[i],
            version[[i]][j],
            "MethodModifications.xsd",
            sep="/"
        )
        xsd <- xml2::read_xml(url)
        .validMsSettingsXsd[[i]][[j]] <- .xmlValidMsSettings(xsd)
    }
}

usethis::use_data(
    .validMsSettingsXsd,
    internal=TRUE,
    overwrite=TRUE
)
