#' Create meth files from xml templates
#'
#' This function calls `XmlMethodChanger.exe` on all given xml files generated
#' with [writeMethodXmls()].
#' It works only on Windows.
#'
#' @param template `character`, path to template `.meth` file.
#' @param xml `character`, vector of path to `.xml` files.
#' @param executable `character`, path to the `XmlMethodChanger.exe`
#' executable.
#' @param verbose `logical`, if `TRUE` a progress bar is shown.
#' @return Nothing. Used for its side effects.
#' @seealso [writeMethodXmls()]
#' @references XmlMethodChanger source code:
#' https://github.com/thermofisherlsms/meth-modifications/
#' @rdname XmlMethodChanger
#' @export
#' @examples
#' \dontrun{
#' runXmlMethodChanger(templateMeth="TMS2IndependentTemplate240Extended.meth",
#'                     modificationXml=list.files(pattern="^method.*\\.xml$"),
#'                     executable="..\\XmlMethodChanger.exe")
#' }
createTngFusionMethFiles <- function(template,
                                     xml=list.files(pattern=".*\\.xml$"),
                                     executable="XmlMethodChanger.exe",
                                     verbose=interactive()) {
    if (.Platform$OS.type != "windows") {
        stop("This function works only on Windows.")
    }
    if (!file.exists(executable)) {
        stop(sQuote(executable), " not found!")
    }
    if (!file.exists(template)) {
        stop(sQuote(template), " not found!")
    }

    if (verbose) {
        pb <- txtProgressBar(0L, length(xml))
    }

    for (i in seq(along=xml)) {
        .xmlMethodChanger(executable, template,
                          .swapFileExt(xml[i]), xml[i])
        if (verbose) {
            setTxtProgressBar(pb, i)
        }
    }

    if (verbose) {
        close(pb)
    }
}

#' @rdname XmlMethodChanger
#' @export
runXmlMethodChanger <- createTngFusionMethFiles

.xmlMethodChanger <- function(exe, input, output, modification) {
    args <- paste0("--", c("input", "output", "modification"), "=",
                   shQuote(c(input, output, modification)))
    system2(normalizePath(exe), args=args)
}

#' Create .csv/.txt files from raw files
#'
#' This function calls `ScanHeadsman.exe` on a given directory containing
#' `.raw` files. `ScanHeadsman.exe` extracts the method and scan header data
#' into `.experiments.csv` and `.txt` files, respectively.
#' It works only on Windows.
#'
#' @param executable `character`, path to the `ScanHeadsman.exe` executable.
#' @param path `character`, path to the directory containing the `.raw` files.
#' @return Nothing. Used for its side effects.
#' @seealso [runXmlMethodChanger()]
#' @rdname ScanHeadsman
#' @references ScanHeadsman source code:
#' https://bitbucket.org/caetera/scanheadsman
#' @export
#' @examples
#' \dontrun{
#' runScanHeadsman("raw", executable="..\\ScanHeadsman.exe")
#' }
runScanHeadsman <- function(path=".", executable="ScanHeadsman.exe") {
    if (.Platform$OS.type != "windows") {
        stop("This function works only on Windows.")
    }
    if (!file.exists(executable)) {
        stop(sQuote(executable), " not found!")
    }
    .scanHeadsman(exe=executable, path=path)
}

.scanHeadsman <- function(exe, path=".") {
    args <- c(paste0("--", c("noMS", "methods:CSV")), path)
    system2(normalizePath(exe), args=args)
}
