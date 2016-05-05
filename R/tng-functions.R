#' @export
createTngFusionMethFiles <- function(templateMeth,
                                     modificationXml=list.files(pattern=".*\\.xml$"),
                                     executable="XmlMethodChanger.exe",
                                     verbose=interactive()) {
  if (.Platform$OS.type != "windows") {
    stop("This function works only on Windows.")
  }
  if (!file.exists(executable)) {
    stop(sQuote(executable), " not found!")
  }
  if (!file.exists(templateMeth)) {
    stop(sQuote(templateMeth), " not found!")
  }

  if (verbose) {
    pb <- txtProgressBar(0, length(modificationXml))
  }

  for (i in seq(along=modificationXml)) {
    .xmlMethodChanger(executable, templateMeth,
                      .swapExtension(modificationXml[i]),
                      modificationXml[i])
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

