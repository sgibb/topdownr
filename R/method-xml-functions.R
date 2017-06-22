#' Create Mass Label
#'
#' Identifying the experiments by the running time/order is complicated.
#' Sometimes the instrument records a new run with the same settings which moves
#' all indicies. Same is true for the start times.
#'
#' We modify the target mass slightly by rounding to one decimal place and use
#' the second to fourth (default) to encode the id. Subsequently it is possible
#' to find the results by their individual mass.
#'
#' @param x double, original mass
#' @param id double, run id
#' @param divisor double, divisor (determines which decimal place
#' @return double, mass label (id encoded in the second to fourth decimal place)
#' @noRd
#' @example
#' library("topdown")
#' topdown:::.massLabel(c(750, 1000), c(1, 100))
.massLabel <- function(x, id, divisor=10000L) {
  if (any(log10(divisor) <= log10(id) + 1L)) {
    stop(sQuote("divisor"), " has to be at least two digits more than ",
         sQuote("id"))
  }
  round(x, digits=1L) + id/divisor
}

#' Create settings for MS2 experiments
#' @param ms2Settings list of MS2 settings, e.g. by defaultMs2Settings()
#' @param groupBy character, colnames used for grouping
#' @param replications integer, how often replicate the settings
#' @parma randomise logical, randomise experiments
#' @return data.frame
#' @noRd
.ms2Experiments <- function(ms2Settings, groupBy=c("replication", "ETDReactionTime"),
                            replications=2L, randomise=TRUE) {
  ms2Settings <- modifyList(ms2Settings, list(replication=1L:replications))
  experiments <- expand.grid(ms2Settings, stringsAsFactors=FALSE)
  experiments <- .replaceZeroETDReactionTime(experiments)
  if (length(groupBy)) {
    experiments <- .groupExperimentsBy(experiments, groupBy)
  } else {
    experiments <- list(experiments)
  }
  if (randomise) {
    experiments <- lapply(experiments, .resample)
  }
  experiments
}

#' Change ETD settings for MS2 experiments
#' @param x data.frame, generated from .ms2Experiments
#' @return data.frame
#' @seealso https://github.com/sgibb/topdown/issues/9
#' @noRd
.replaceZeroETDReactionTime <- function(x) {
  col <- grep("ETDReactionTime", colnames(x))
  if (length(col)) {
    isZero <- x[, col] == 0L
    activationType <- toupper(gsub("^ET", "", x$ETDSupplementalActivation[isZero]))
    x$ActivationType[isZero] <- activationType
    x$CollisionEnergy[isZero] <- x$ETDSupplementalActivationEnergy[isZero]
    x[isZero, grepl("^ETD", colnames(x))] <- NA
  }
  x
}

#' Split MS2 experiment data.frame
#' @param x data.frame, from .ms2Experiments
#' @param cols character, colnames used to split
#' @return list
#' @noRd
.groupExperimentsBy <- function(x, cols=c("replication")) {
  if (length(cols) > 1L) {
    ## `interaction` doesn't handle NA values, so use `paste` instead
    f <- do.call(paste, c(x[, cols], sep=":"))
  } else {
    f <- x[, cols]
  }
  split(x, f)
}

#' Calculate Start/End time of each scan
#' @param nMs2 integer, number of MS2 scans
#' @param nMs2perMs1 integer, number of MS2 scans for each MS1 scan
#' @param duration double, duration of scan in minutes
#' @param gap double, pause between two scans in minutes
#' @return data.frame with columns type (type of scan), start and end times in
#' minutes
.startEndTime <- function(nMs2, nMs2perMs1=9L, duration=0.5, gap=0.01) {
  nMs1 <- ceiling(nMs2 * 1L/nMs2perMs1)
  n <- nMs2 + nMs1

  if (n > 300) {
    warning("More than 300 experiments might cause the MS device ",
            "to become unresponsive. Choose other ", sQuote("groupBy"),
            " parameters to reduce number of experiments per file.")
  }

  start <- seq(gap, by=duration, length.out=n)
  end <- seq(duration, by=duration, length.out=n)

  type <- rep_len(rep.int(c("MS1", "MS2"), times=c(1L, nMs2perMs1)), n)
  data.frame(Type=type, StartTimeMin=start, EndTimeMin=end,
             stringsAsFactors=FALSE)
}

#' Resample rows (Experiments) in a data.frame
#' @param x data.frame
#' @param fun function to apply on the number of rows, e.g. use `seq` to do
#' nothing
#' @return reordered data.frame
#' @noRd
.resample <- function(x, fun=sample) {
  fun <- match.fun(fun)
  x[fun(nrow(x)),]
}

#' Valid XML MS1 tags
#' @noRd
.validMs1Tags <- function() {
  c("FirstMass", "LastMass", "Microscans", "MaxITTimeInMS", "AgcTarget")
}

#' Valid XML MS2 tags
#' @noRd
.validMs2Tags <- function() {
  c("ActivationType", "IsolationWindow", "EnableMultiplexIons",
    "EnableMSXIds", "MaxNoOfMultiplexIons", "OrbitrapResolution",
    "AgcTarget", "MinAgcTarget", "MaxITTimeInMS", "Microscans",
    "ETDReactionTime", "ETDReagentTarget",
    "MaximumETDReagentInjectionTime", "UseInternalCalibratedETD",
    "ETDSupplementalActivationEnergy", "ETDSupplementalActivation")
}

#' XML writing functions
#' @param file filename
#' @param encoding file encoding
#' @noRd
.xmlHeader <- function(file, encoding="utf-8") {
  cat0("<?xml version=\"1.0\" encoding=\"", encoding, "\"?>\n", file=file)
}

#' @param name character, tag name
#' @param value tag value
#' @param attrs named vector, names == attribute names
#' @param close logical, should the tag be closed
#' @param indention integer, number of spaces used for indention
#' @param file filename
#' noRd
.xmlTag <- function(name, value=character(), attrs=character(), close=TRUE,
                    indention=0L, file) {
  indention <- paste0(rep(" ", times=indention))

  if (length(attrs)) {
    attrs <- paste0(" ", names(attrs), "=\"", attrs, "\"")
  }
  if (length(value) && close) {
    cat0(indention, "<", name, attrs, ">", value, "</", name, ">\n", file=file)
  } else if (!length(value) && close) {
    cat0(indention, "<", name, attrs, "/>\n", file=file)
  } else {
    cat0(indention, "<", name, attrs, ">", value, "\n", file=file)
  }
}

#' @param name character, tag name
#' @param indention integer, number of spaces used for indention
#' @param file filename
#' noRd
.xmlTagClose <- function(name, indention=0L, file) {
  cat0(paste0(rep(" ", times=indention)), "</", name, ">\n", file=file)
}

#' @param x list, named
#' @param indention integer, number of spaces used for indention
#' @param na.rm should NAs be removed?
#' @param file filename
#' noRd
.xmlListToTags <- function(x, indention=0L, na.rm=TRUE, file) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  invisible(mapply(.xmlTag, name=names(x), value=x,
                   MoreArgs=list(indention=indention, file=file),
                   SIMPLIFY=FALSE, USE.NAMES=FALSE))
}

#' Write method.xml file.
#' @param ms1 list, ms1 settings
#' @param ms2 data.frame, ms2 experiment settings
#' @param times data.frame, experiment times
#' @param massList matrix, 2 columns (mz, z)
#' @param massLabeling logical, should mz values modified to used as labels?
#' @param file filename
#' @param encoding file encoding
#' noRd
.writeMethodXml <- function(ms1, ms2, times, massList, massLabeling=TRUE,
                            file, encoding="utf-8") {
  ## stop if file isn't writeable
  if (file.exists(file) && file.access(file, 2) != 0) {
    stop("No permissions to write into ", sQuote(file), "!")
  }

  ## file handle
  f <- file(file, open="wt", encoding=encoding)
  on.exit(close(f))

  ## header
  .xmlHeader(file=f, encoding=encoding)

  ## MethodModifications frame
  .xmlTag("MethodModifications", attrs=c(Version=1,
                                         Model="OrbitrapFusion",
                                         Family="Calcium",
                                         Type="SL"), close=FALSE, file=f)

  ## First MS1 scan
  .xmlFullMsScan(ms1, file=f)

  n <- nrow(times)

  ## Copy experiments
  for (i in 3L:n) {
    .xmlCopyAndAppendExperiment(order=i - 1L,
                                src=as.integer(times$Type[i] == "MS2"),
                                file=f)
  }

  ## Start/EndTime
  for (i in 1L:n) {
    .xmlStartEndTime(order=n + i - 1L,
                     times=c(times[i, c("StartTimeMin", "EndTimeMin")]),
                     file=f)
  }

  ## mass labeling
  nms2 <- nrow(ms2)
  mz <- replicate(nms2, massList[, 1L])
  if (massLabeling) {
    mz <- .massLabel(mz, id=rep(1L:nms2, each=nrow(mz)))
  }

  ## TMSn scans
  for (i in 1L:nms2) {
    .xmlTMSnScan(ms2[i, ], mz=mz[, i], z=massList[, 2L],
                 order=2L * n + i - 1L, idx=i, file=f)
  }

  .xmlTagClose("MethodModifications", file=f)
}

#' Full MS Scan tag
#' @param x MS1 settings
#' @param file filename
#' @noRd
.xmlFullMsScan <- function(x, file) {
  .xmlTag("Modification", attrs=c(Order=1), close=FALSE, file=file)
  .xmlTag("Experiment", attrs=c(ExperimentIndex=0), close=FALSE, indention=2L,
          file=file)
  .xmlTag("FullMSScan", close=FALSE, indention=4L, file=file)
  .xmlListToTags(x, indention=6L, file=file)
  .xmlTagClose("FullMSScan", indention=4L, file=file)
  .xmlTagClose("Experiment", indention=2L, file=file)
  .xmlTagClose("Modification", file=file)
}

#' TMSn Scan tag
#' @param x data.frame, with one row, MS2 settings
#' @param mz double
#' @param z double
#' @param order integer, modification order
#' @param idx integer, experiment index
#' @param file filename
#' @noRd
.xmlTMSnScan <- function(x, mz, z, order, idx, file) {
  .xmlTag("Modification", attrs=c(Order=order), close=FALSE, file=file)
  .xmlTag("Experiment", attrs=c(ExperimentIndex=idx), close=FALSE,
          indention=2L, file=file)
  .xmlTag("TMSnScan", close=FALSE, indention=4L, file=file)
  .xmlListToTags(c(x[, !grepl("CollisionEnergy", colnames(x))]), indention=6L,
                 file=file)
  .xmlMassList(mz=mz, z=z, energy=x$CollisionEnergy, type=x$ActivationType,
               indention=6L, file=file)
  .xmlTagClose("TMSnScan", indention=4L, file=file)
  .xmlTagClose("Experiment", indention=2L, file=file)
  .xmlTagClose("Modification", file=file)
}

#' Copy and AppendExperiment tag
#' @param order integer, number of experiments
#' @param src integer, source experiment index
#' @param file filename
#' @noRd
.xmlCopyAndAppendExperiment <- function(order, src, file) {
  .xmlTag("Modification", attrs=c(Order=order), close=FALSE, file=file)
  .xmlTag("CopyAndAppendExperiment", attrs=c(SourceExperimentIndex=src),
          indention=2L, file=file)
  .xmlTagClose("Modification", file=file)
}

#' Start/EndTimeMin tag
#' @param order integer, number of experiments
#' @param times double, named vector with times
#' @param file filename
#' @noRd
.xmlStartEndTime <- function(order, times, file) {
  stopifnot(length(times) == 2)
  .xmlTag("Modification", attrs=c(Order=order), close=FALSE, file=file)
  .xmlListToTags(setNames(times, c("StartTimeMin", "EndTimeMin")),
                 indention=2, file=file)
  .xmlTagClose("Modification", file=file)
}

#' MassList tag
#' @param mz double
#' @param z double
#' @param energy double, length=1
#' @param type activation type
#' @param indention integer, number of spaces used for indention
#' @param file filename
.xmlMassList <- function(mz, z, energy=0L, type=c("ETD", "HCD", "CID"),
                         indention=0L, file) {
  type <- match.arg(type)
  if (type != "ETD") {
    .xmlTag("MassList", attrs=setNames("true", paste0("CollisionEnergy", type)),
            close=FALSE, indention=indention, file=file)
  } else {
    .xmlTag("MassList", close=FALSE, indention=indention, file=file)
  }
  for (i in seq(along=mz)) {
    .xmlMassListRecord(mz=mz[i], z=z[i], energy=energy, type=type,
                       indention=indention + 2L, file=file)
  }
  .xmlTagClose("MassList", indention=indention, file=file)
}

#' MassListRecord tag
#' @param mz double, length=1
#' @param z double, length=1
#' @param energy double, length=1
#' @param type activation type
#' @param indention integer, number of spaces used for indention
#' @param file filename
.xmlMassListRecord <- function(mz, z, energy=0L, type=c("ETD", "HCD", "CID"),
                               indention=0L, file) {
  type <- match.arg(type)
  .xmlTag("MassListRecord", close=FALSE, indention=indention, file=file)
  .xmlTag("MOverZ", value=mz, indention=indention + 2L, file=file)
  ## z == 10 is the maximal allowed charge state
  .xmlTag("Z", value=min(z, 10L), indention=indention + 2L, file=file)
  if (type != "ETD") {
    .xmlTag(paste0("CollisionEnergy", type), value=energy,
            indention=indention + 2L, file=file)
  }
  .xmlTagClose("MassListRecord", indention=indention, file=file)
}



#'@export
writeMethodXmls <- function(ms1Settings, ms2Settings, replications=2,
                            groupBy=c("replication",
                                      "ETDReactionTime"),
                            massList=NULL,
                            nMs2perMs1=10, duration=0.5,
                            randomise=TRUE, pattern="method_%s.xml") {

  ## TODO: test for valid ms1 and ms2 tags before creating files

  mods <- .modifications(ms1Settings=ms1Settings,
                         ms2Settings=ms2Settings,
                         replications=replications,
                         groupBy=groupBy, massList=massList,
                         nMs2perMs1=nMs2perMs1, duration=duration,
                         randomise=randomise)

  files <- sprintf(pattern, names(mods$ms2))

  for (i in seq(along=mods$ms2)) {
    cat(saveXML(.toMethodModificationXml(ms1=mods$ms1, ms2=mods$ms2[[i]],
                                         times=mods$times,
                                         massList=mods$massList),
                encoding="utf8"), file=files[i])
  }
}
