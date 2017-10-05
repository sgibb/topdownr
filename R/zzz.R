.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste(
            "\nThis is topdownr version", packageVersion("topdownr"), "\n",
            " Visit https://sgibb.github.io/topdownr/ to get started.\n"))
}
