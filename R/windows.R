#' Windows specific functions.
#'
#' The functions `runXmlMethodChanger` and `runScanHeadsman` call
#' `XmlMethodChanger.exe` and `ScanHeadsman.exe` with the correspond arguments.
#' The only work on Windows (maybe on Linux + wine as well but that was never
#' tested).
#'
#' @details
#'
#' `runXmlMethodChanger` applies ‘XmlMethodChanger.exe’ on all given XML files
#' generated with [writeMethodXmls()] to create `.meth` files from a template.
#'
#' @param template `character`, path to template `.meth` file.
#' @param xml `character`, vector of path to `.xml` files.
#' @param executable `character`, path to the `XmlMethodChanger.exe` or
#' `ScanHeadsman.exe` executable.
#' @param verbose `logical`, if `TRUE` a progress bar is shown.
#' @return Nothing. Used for its side effects.
#' @seealso [writeMethodXmls()]
#' @references XmlMethodChanger source code:
#' https://github.com/thermofisherlsms/meth-modifications/
#' @rdname windows-specific-functions
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

#' @rdname windows-specific-functions
#' @export
runXmlMethodChanger <- createTngFusionMethFiles

.xmlMethodChanger <- function(exe, input, output, modification) {
    args <- paste0("--", c("input", "output", "modification"), "=",
                   shQuote(c(input, output, modification)))
    system2(normalizePath(exe), args=args)
}

#' @rdname windows-specific-functions
#'
#' @details
#'
#' `runScanHeadsman` calls `ScanHeadsman.exe` on a given directory containing
#' `.raw` files. `ScanHeadsman.exe` extracts the method and scan header data
#' into `.experiments.csv` and `.txt` files, respectively.
#'
## @param executable `character`, path to the `ScanHeadsman.exe` executable.
#' @param path `character`, path to the directory containing the `.raw` files.
## @return Nothing. Used for its side effects.
#' @export
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
