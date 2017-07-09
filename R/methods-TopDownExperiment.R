#' @param x TopDownExperiment
#' @noRd
setMethod("[", "TopDownExperiment", function(x, i, j="missing", ...,
                                             drop="missing") {

  if (missing(i)) {
    i <- seq_along(x)
  }
  if (!(is.logical(i) | is.numeric(i) | is.character(i))) {
    stop("subsetting works only with numeric, logical or character")
  }
  if (length(x) < length(i)) {
    stop("subscript out of bounds")
  }
  if (is.character(i)) {
    i <- unique(match(i, featureNames(x)))
  }
  if (anyNA(i)) {
    stop("subscript out of bounds")
  }
  if (is.numeric(i)) {
    if (max(i) > length(x) | min(i) < 1) {
      stop("subscript out of bounds")
    }
  }
  d0 <- dim(x)
  fn <- featureNames(x)[i]
  x@assignmentTable <- copy(assignmentTable(x)[SpectrumId %in% fn,])

  ## overwrite MSnbase log msg (and duplicated log by TDE subsetting if j is
  ## given)
  processing <- x@processingData@processing
  x <- callNextMethod(x=x, i=i, j=j, ..., drop=drop)

  if (!missing(j) && length(j) && j[1L] != "missing") {
    if (is.character(j)) {
      x <- .filterFragmentIonOrType(x, j)
    } else if (is.numeric(j)) {
      x <- .filterFragmentPos(x, j)
    } else if (is.logical(j) && length(j) == nrow(fragmentTable(x))) {
      x <- .filterFragmentId(x, .fragmentId(x)[j])
    } else {
      stop("subscript out of bounds")
    }
  }

  d1 <- dim(x)
  x@processingData@processing <- processing
  x <- .logmsg(x, paste0("Subsetted [", d0[1L], ";", d0[2L], "] to [",
                                        d1[1L], ";", d1[2L], "]."))

  x
})

#' aggregate TopDownExperiment objects
#' @param x TopDownExperiment
#' @param by character, fvarLabels used for aggregation
#' @param verbose logical, verbose output?
#' @return TopDownExperiment
#' @noRd
setMethod("aggregate", "TopDownExperiment", function(x,
  by=c("Mz", "AGCTarget", "ETDReagentTarget", "ETDActivation", "CIDActivation",
       "HCDActivation", "SupplementalActivationCE", "SupplementalActivation"),
  ..., verbose=interactive()) {

  if (!(all(by %in% fvarLabels(x)))) {
    stop(by[!by %in% fvarLabels(x)], " is/are not present in 'fvarLabels(x)'.")
  }

  fd <- as(featureData(x), "data.frame")
  group <- .groupByLabels(fd, by)

  .msg(verbose, "Aggregating featureData.")

  fd <- .aggregateDataFrame(fd, f=group,
                            ignoreCols=c("Scan", "spectrum", "ConditionId", "File"))
  mz <- split(mz(x), group)
  int <- split(intensity(x), group)
  rtm <- split(rtime(x), group)
  fns <- split(featureNames(x), group)
  newFns <- MSnbase:::formatFileSpectrumNames(0L, seq_along(mz),
                                              nSpectra=length(mz), nFiles=0L)

  if (verbose) {
    message("Aggregating spectra.")
    pb <- txtProgressBar(min=0L, max=length(fns), style=3L)
  }

  atab <- copy(x@assignmentTable)
  newAssay <- new.env(parent=emptyenv())

  for (i in seq(along=fns)) {
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    assign(newFns[i],
           .aggregateSpectra(mz[[i]], int[[i]], rt=rtm[[i]],
                             atab[SpectrumId %in% fns[[i]],],
                             fragmentTable(x), acquisitionNum=i, fromFile=1L),
           newAssay)
    atab[SpectrumId %in% fns[[i]], SpectrumId := newFns[i]]
  }

  if (verbose) {
    close(pb)
  }

  atab[, MzId := NULL]
  setorder(atab, SpectrumId, FragmentId)
  atab <- unique(atab)
  atab <- .updateAssignmentTableMzId(atab)

  lockEnvironment(newAssay, bindings=TRUE)

  rownames(fd) <- newFns
  fdata <- new("AnnotatedDataFrame", data=fd)

  pd <- data.frame(sampleNames="aggregated")
  rownames(pd) <- pd$sampleNames
  pdata <- new("NAnnotatedDataFrame", data=pd)

  .msg(verbose, "Creating aggregated object.")

  td <- new("TopDownExperiment",
            assayData=newAssay,
            featureData=fdata,
            phenoData=pdata,
            experimentData=experimentData(x),
            processingData=processingData(x),
            sequence=x@sequence,
            assignmentTable=atab,
            fragmentTable=fragmentTable(x))

  d <- c(dim(x), dim(td))
  td <- .logmsg(td, paste0("Aggregated [", d[1L], ";", d[2L], "] to [",
                                           d[3L], ";", d[4L], "]."))

  if (validObject(td)) {
    td
  }
})

#' Accessor for assignmentTable, not exported yet
#' @param object TopDownExperiment
#' @return data.table
#' @noRd
setMethod("assignmentTable", "TopDownExperiment", function(object) {
  object@assignmentTable
})

#' Dim (number of spectra, number of spectra-fragment-assignments)
#' @param object TopDownExperiment
#' @return numeric
#' @noRd
setMethod("dim", "TopDownExperiment", function(x) {
  as.integer(c(length(x), nrow(assignmentTable(x))))
})

#' Accessor for fragmentTable, not exported yet
#' @param object TopDownExperiment
#' @return data.table
#' @noRd
setMethod("fragmentTable", "TopDownExperiment", function(object) {
  object@fragmentTable
})

#' @param object TopDownExperiment
#' @noRd
setMethod("show", "TopDownExperiment", function(object) {
  cat("TopDown Experiment data (\"", class(object), "\").\n", sep="")

  if (nchar(object@sequence)) {
    cat("- - - Protein data - - -\n")
    prefix <- sprintf("Amino acid sequence (%d):", nchar(object@sequence))
    cat(prefix, .snippet(object@sequence, getOption("width") - nchar(prefix)),
        "\n")
  }

  if (nrow(fragmentTable(object))) {
    cat("- - - Fragment data - - -\n")
    cat("Number of theoretical fragments:", nrow(fragmentTable(object)), "\n")
    fragments <- sort(unique(.fragmentTypes(object)))
    prefix <- sprintf("Theoretical fragment types (%d):", length(fragments))
    cat(prefix, .snippet(paste0(fragments, collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (nrow(assignmentTable(object))) {
    cat("- - - Spectra data - - -\n")
    cat("Number of Spectra:", length(object), "\n")
    cat("Number of Spectra-Fragment-Assignments:",
        nrow(assignmentTable(object)), "\n")
    fragments <- sort(unique(.fragmentTypes(object)[
      .fragmentId(object) %in% assignmentTable(object)$FragmentId]))
    prefix <- sprintf("Assigned fragment types (%d):", length(fragments))
    cat(prefix, .snippet(paste0(fragments, collapse=", "),
                         getOption("width") - nchar(prefix)), "\n")
  }

  if (length(processingData(object)@processing)) {
    cat("- - - Processing information - - -\n")
    cat(paste0(processingData(object)@processing, collapse="\n"), "\n")
  }

  invisible(NULL)
})

