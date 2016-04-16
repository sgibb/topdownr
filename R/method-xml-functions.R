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
  xml$addNode("Modification", attrs=c(Order=order), ...)
  xml
}

.xmlExperiment <- function(xml, index, ...) {
  xml$addNode("Experiment", attrs=c(ExperimentIndex=index), ...)
  xml
}

.xmlCopyAndAppendExperiment <- function(xml, order, source) {
  xml <- .xmlModification(xml, order=order, close=FALSE)
  xml$addTag("CopyAndAppendExperiment", attrs=c(SourceExperimentIndex=source))
  xml$closeTag()
  xml
}

.xmlFullMSScan <- function(xml, index, settings=defaultMs1Settings()) {
  xml <- .xmlExperiment(xml, index=index, close=FALSE)
  xml$addNode("FullMSScan",
              .children=mapply(xml$addTag,
                               name=names(settings), value=settings,
                               SIMPLIFY=FALSE, USE.NAMES=FALSE))
  xml$closeTag()
  xml
}

writeMethodXml <- function(file=NULL, ms1Settings=defaultMs1Settings(), ...) {
  if (!is.null(file) && file.exists(file)) {
    stop("File already exists!")
  }

  xml <- .xmlMethodModification()
  xml <- .xmlFullMSScan(xml, index=0, settings=ms1Settings)

  cat(saveXML(xml, encoding="utf8"), file=file)
}
