#' Create meth files from xml templates
#'
#' This function calls \sQuote{XmlMethodChanger.exe} on all given xml
#' files generated with \code{\link{writeMethodXmls}}.
#' It works only on Windows.
#'
#' @param template \code{character}, path to template meth file.
#' @param xml \code{character}, vector of path to xml files.
#' @param executable \code{character}, path to the
#'  \sQuote{XmlMethodChanger.exe} executable.
#' @param verbose \code{logical}, if \code{TRUE} a progress bar is shown.
#' @return Nothing. Used for its side effects.
#' @seealso \code{\link{writeMethodXmls}}
#' @references XmlMethodChanger source code:
#' \url{https://github.com/thermofisherlsms/meth-modifications}
#' @rdname XmlMethodChanger
#' @export
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
#' This function calls \sQuote{ScanHeadsman.exe} on a given directory containing
#' RAW files.
#' It works only on Windows.
#'
#' @param executable \code{character}, path to the
#'  \sQuote{ScanHeadsman.exe} executable.
#' @param path \code{character}, path to the directory containing the raw files.
#' @return Nothing. Used for its side effects.
#' @seealso \code{\link{XmlMethodChanger}}
#' @rdname ScanHeadsman
#' @references ScanHeadsman source code:
#' \url{https://bitbucket.org/caetera/scanheadsman}
#' @export
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

