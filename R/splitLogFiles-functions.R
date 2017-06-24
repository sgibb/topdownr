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
