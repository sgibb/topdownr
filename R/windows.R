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
                      .swapExtension(xml[i]), xml[i])
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
runScanHeadsman <- function(executable="ScanHeadsman.exe",
                            path=".") {
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

.swapExtension <- function(x, ext="meth") {
  paste(file_path_sans_ext(x), ext, sep=".")
}

#' Split Log files.
#'
#' This function watches a Xtract logfile and split it into individual spectra
#' files.
#'
#' @param file file to watch
#' @param idletime time to wait for file changes
#' @param sleep time to sleep before next check
#' @param verbose verbose output?
#' @export
splitLogFile <- function(file, idletime=20L, sleep=0.1,
                         verbose=interactive()) {
  startTime <- proc.time()[3L]
  mTime <- 0L

  while (proc.time()[3L] - startTime < idletime) {
    mTimeCur <- file.mtime(file)
    if (file.exists(file) && mTimeCur > mTime && .revfgrepl("</Xtract>", file)) {
      nScan <- .scanNumber(file)
      newFile <- paste0(file_path_sans_ext(file), nScan, ".xml")
      if (verbose) {
        message("Renaming ", sQuote(file), " to " , sQuote(newFile), ".")
      }
      file.rename(file, newFile)
      mTime <- mTimeCur
      startTime <- proc.time()[3L]
    }
    Sys.sleep(sleep)
  }
  invisible(NULL)
}

#' reads just the last characters of a file needed to match against
#' @param pattern pattern to look for, no regex supported
#' @param file file name
#' @noRd
.revfgrepl <- function(pattern, file) {
  size <- file.size(file)
  n <- nchar(pattern)

  if (size < n) {
    return(FALSE)
  }

  f <- file(file, "rb")
  on.exit(close(f))
  seek(f, where=size-n)
  pattern == readChar(f, nchars=n)
}

#' fetch the scan number
#' @param file file name
#' @param nchars read the first nchars characters
#' @noRd
.scanNumber <- function(file, nchars=1024L) {
  f <- file(file, "rb")
  on.exit(close(f))
  ch <- readChar(f, nchars=nchars)
  rx <- regexpr("(?<=Scan=\")\\d+", ch, perl=TRUE)
  as.double(regmatches(ch, rx))
}
