if (is.null(getGeneric("assayData"))) {
  setGeneric("assayData", function(object) standardGeneric("assayData"))
}
if (is.null(getGeneric("rowViews"))) {
  setGeneric("rowViews", function(object, ...) standardGeneric("rowViews"))
}
if (is.null(getGeneric("fragmentViews"))) {
  setGeneric("fragmentViews", function(object, ...) standardGeneric("fragmentViews"))
}
if (is.null(getGeneric("colData"))) {
  setGeneric("colData", function(object, ...) standardGeneric("colData"))
}
if (is.null(getGeneric("conditionData"))) {
  setGeneric("conditionData", function(object, ...) standardGeneric("conditionData"))
}
