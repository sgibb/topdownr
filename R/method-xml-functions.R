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
#' @param replications integer, how often replicate the settings
#' @return data.frame
#' @noRd
.ms2Experiments <- function(ms2Settings, replications=2L) {
  ms2Settings <- modifyList(ms2Settings, list(replication=1L:replications))
  expand.grid(ms2Settings, stringsAsFactors=FALSE)
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
#' @return list with multiple (splitted) data.frame(s)
#' @noRd
.groupExperimentsBy <- function(x, cols=c("replication")) {
  if (length(cols) > 1L) {
    f <- interaction(as.list(x[, cols]))
  } else {
    f <- x[, cols]
  }
  split(x, f)
}

.startEndTime <- function(nMs2, nMs2perMs1=9, duration=0.5, gap=0.01) {
  nMs1 <- ceiling(nMs2 * 1/nMs2perMs1)
  n <- nMs2 + nMs1

  start <- seq(gap, by=duration, length.out=n)
  end <- seq(duration, by=duration, length.out=n)

  type <- rep_len(rep.int(c("MS1", "MS2"), times=c(1L, nMs2perMs1)), n)
  data.frame(type=type, StartTimeMin=start, EndTimeMin=end)
}

## could use `seq` too for ordered output
.resample <- function(x, fun=sample) {
  fun <- match.fun(fun)
  x[fun(nrow(x)),]
}

.modifications <- function(ms1Settings, ms2Settings, replications=2,
                           groupBy=c("replication",
                                     "ETDReactionTime"),
                           massList=NULL,
                           nMs2perMs1=10, duration=0.5,
                           randomise=TRUE) {

  ms2Experiments <- .ms2Experiments(ms2Settings, replications)
  ms2Experiments <- .replaceZeroETDReactionTime(ms2Experiments)
  ms2Experiments <- .groupExperimentsBy(ms2Experiments, groupBy)
  times <- .startEndTime(max(sapply(ms2Experiments, nrow)), nMs2perMs1)
  if (nrow(times) > 300) {
    warning("More than 300 experiments might cause the MS device ",
            "to become unresponsive. Choose other ", sQuote("groupBy"),
            " parameters to reduce number of experiments per file.")
  }
  if (randomise) {
    ms2Experiments <- lapply(ms2Experiments, .resample)
  }
  list(ms1=ms1Settings,
       ms2=ms2Experiments,
       times=times,
       massList=massList)
}

.toMethodModificationXml <- function(ms1, ms2, times, massList) {
  xml <- .xmlMethodModification()
  xml <- .xmlScanTemplates(xml, ms1, ms2, massList)

  nr <- nrow(times)
  times$id <- (1:nr) - 1
  times$ms2idx[times$type == "MS2"] <- 1:nrow(ms2)

  # 1 and 2 are the two template scans
  ## copy templates
  for (i in 3:nr) {
    xml <- .xmlCopyAndAppendExperiment(xml, order=i,
                                       source=as.numeric(times$type[i] == "MS2"))
  }

  ## change start/end times
  for (i in 1:nr) {
    xml <- .xmlModification(xml, order=i + nr, close=FALSE)
    xml <- .xmlExperiment(xml, index=times$id[i], close=FALSE)
      xml <- .xmlListToTags(xml, as.list(times[i, c("StartTimeMin", "EndTimeMin")]))
    xml$closeTag()
    xml$closeTag()
  }

  ## modify templates
  xml <- .xmlModification(xml, order=2 * nr + 1, close=FALSE)
  for (i in 1:nr) {
    if (times$type[i] == "MS2") {
      xml <- .xmlScan(xml, index=times$id[i], type="TMSnScan",
                      settings=ms2[times$ms2idx[i],
                                   colnames(ms2) %in% .validMs2Tags()])
    }
  }

  xml$closeTag()

  xml
}

.xmlCopyAndAppendExperiment <- function(xml, order, source) {
  xml <- .xmlModification(xml, order=order, close=FALSE)
  xml$addTag("CopyAndAppendExperiment", attrs=c(SourceExperimentIndex=source))
  xml$closeTag()
  xml
}

.xmlMethodModification <- function(Model="OrbitrapFusion", Family="Calcium",
                                   Type="SL") {
  xml <- withCallingHandlers(xmlTree("MethodModifications",
                                     attrs=c(Version=1,
                                             Model="OrbitrapFusion",
                                             Family="Calcium",
                                             Type="SL")),
                             warning=function(w) {
                               if (any(grepl("empty XML", w))) {
                                 invokeRestart("muffleWarning")
                               }
                             })
  xml
}

.xmlModification <- function(xml, order, ...) {
  xml$addTag("Modification", attrs=c(Order=order), ...)
  xml
}

.xmlExperiment <- function(xml, index, ...) {
  xml$addTag("Experiment", attrs=c(ExperimentIndex=index), ...)
  xml
}

.xmlCopyAndAppendExperiment <- function(xml, order, source) {
  xml <- .xmlModification(xml, order=order, close=FALSE)
  xml$addTag("CopyAndAppendExperiment", attrs=c(SourceExperimentIndex=source))
  xml$closeTag()
  xml
}

.xmlScan <- function(xml, index, settings, type=c("FullMSScan", "TMSnScan"),
                     close=TRUE) {
  type <- match.arg(type)
  xml <- .xmlExperiment(xml, index=index, close=FALSE)
    xml$addTag(type, close=FALSE)
    xml <- .xmlListToTags(xml, settings)
  if (close) {
    xml$closeTag() # Scan
    xml$closeTag() # Experiment
  }
  xml
}

.xmlListToTags <- function(xml, l) {
  mapply(xml$addTag, name=names(l), value=l, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  xml
}

.xmlMassList <- function(xml, l) {
  xml$addTag("MassList", close=FALSE)
    lapply(l, .xmlMassListRecord, xml=xml)
  xml$closeTag()
  xml
}

.xmlMassListRecord <- function(xml, r) {
  xml$addTag("MassListRecord", close=FALSE)
    for (i in 1:nrow(r)) {
      xml$addTag("MOverZ", r[i, "mass"], sep="")
      ## z == 10 is the maximal allowed charge state
      xml$addTag("Z", min(r[i, "z"], 10L), sep="")
    }
  xml$closeTag()
  xml
}

.xmlScanTemplates <- function(xml, ms1Settings, ms2Settings,
                              massList=list()) {
  xml <- .xmlModification(xml, order=1L, close=FALSE)
    xml <- .xmlScan(xml, index=0L, settings=ms1Settings,
                    type="FullMSScan")
  xml$closeTag()

  xml <- .xmlModification(xml, order=2L, close=FALSE)
    xml <- .xmlScan(xml, index=1L, settings=ms2Settings[, colnames(ms2Settings) %in% .validMs2Tags()],
                    type="TMSnScan", close=!length(massList))
    if (length(massList)) {
      xml <- .xmlMassList(xml, massList)
      xml$closeTag() # TMSnScan
      xml$closeTag() # Experiment
    }
  xml$closeTag()
  xml
}

#'@export
writeMethodXmls <- function(ms1Settings, ms2Settings, replications=2,
                            groupBy=c("replication",
                                      "ETDReactionTime"),
                            massList=NULL,
                            nMs2perMs1=10, duration=0.5,
                            randomise=TRUE, pattern="method_%s.xml") {

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
