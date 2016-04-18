.ms2Experiments <- function(ms2Settings, replications=2L) {
  ms2Settings <- modifyList(ms2Settings, list(replication=1:replications))
  expand.grid(ms2Settings, stringsAsFactors=FALSE)
}

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
  data.frame(type=type, StartTime=start, EndTime=end)
}

.replaceZeroETDReactionTime <- function(x) {
  col <- grep("ETDReactionTime", colnames(x))
  if (length(col)) {
    isZero <- x[, col] == 0
    x$ActivationType[isZero] <- toupper(
                                  substr(x$ETDSupplementalActivation[isZero],
                                         3, 5))
    x$CIDCollisionEnergy[isZero] <- x$ETDSupplementalActivationEnergy[isZero]
    x[isZero, grepl("^ETD", colnames(x))] <- NA
  }
  x
}

## could use `seq` to for ordered output
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

  # 1 and 2 are the two template scans
  order <- (1L:nrow(times)) + 2L
  times$id <- 1:nrow(times) - 1L
  times$sources <- ifelse(times$type == "MS1", 0L, 1L)
  times$ms2idx[times$type == "MS2"] <- 1:nrow(ms2)

  ## copy templates
  for (i in seq(along=order)) {
    xml <- .xmlCopyAndAppendExperiment(xml,
                                       order=order[i],
                                       source=times$sources[i])
  }

  order <- order[length(order)] + 1L
  xml <- .xmlModification(xml, order=order, close=FALSE)

  for (i in 1:nrow(times)) {
    if (times$type[i] == "MS1") {
      xml <- .xmlExperiment(xml, index=times$id[i], close=FALSE)
      xml <- .xmlListToTags(xml, as.list(times[i, c("StartTime", "EndTime")]))
      xml$closeTag()
    } else {
      xml <- .xmlScan(xml, index=times$id[i], type="TMSnScan",
                      settings=c(as.list(times[i, c("StartTime", "EndTime")]),
                                 ms2[times$ms2idx[i],]))
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
    xml <- .xmlScan(xml, index=1L, settings=ms2Settings,
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
