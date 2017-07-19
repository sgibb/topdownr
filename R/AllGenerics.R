if (is.null(getGeneric("assayData"))) {
  setGeneric("assayData",
             function(object)standardGeneric("assayData"))
}
if (is.null(getGeneric("rowViews"))) {
  setGeneric("rowViews",
             function(object, ...)standardGeneric("rowViews"))
}
if (is.null(getGeneric("fragmentData"))) {
  setGeneric("fragmentData",
             function(object, ...)standardGeneric("fragmentData"))
}
if (is.null(getGeneric("colData"))) {
  setGeneric("colData",
             function(object, ...)standardGeneric("colData"))
}
if (is.null(getGeneric("colData<-"))) {
  setGeneric("colData<-",
             function(object, ..., value)standardGeneric("colData<-"))
}
if (is.null(getGeneric("conditionData"))) {
  setGeneric("conditionData",
             function(object, ...)standardGeneric("conditionData"))
}
if (is.null(getGeneric("conditionData<-"))) {
  setGeneric("conditionData<-",
             function(object, ..., value)standardGeneric("conditionData<-"))
}
if (is.null(getGeneric("filterIntensity"))) {
  setGeneric("filterIntensity",
             function(object, minIntensity)standardGeneric("filterIntensity"))
}
