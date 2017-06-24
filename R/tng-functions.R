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
                      .swapExtension(xml[i]), xml[i])
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

  if (verbose) {
    close(pb)
  }
}

.xmlMethodChanger <- function(exe, input, output, modification) {
  args <- paste0("--", c("input", "output", "modification"), "=",
                 shQuote(c(input, output, modification)))
  system2(normalizePath(exe), args=args)
}

.swapExtension <- function(x, ext="meth") {
  paste(file_path_sans_ext(x), ext, sep=".")
}

