if (is.null(getGeneric("assayData"))) {
    setGeneric("assayData",
        function(object)standardGeneric("assayData"))
}
if (is.null(getGeneric("rowViews"))) {
    setGeneric("rowViews",
        function(object, ...)standardGeneric("rowViews"))
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
if (is.null(getGeneric("conditionNames"))) {
    setGeneric("conditionNames",
        function(object)standardGeneric("conditionNames"))
}
if (is.null(getGeneric("bestConditions"))) {
    setGeneric("bestConditions",
        function(object, ...)standardGeneric("bestConditions"))
}
if (is.null(getGeneric("filterCv"))) {
    setGeneric("filterCv",
        function(object, ...)standardGeneric("filterCv"))
}
if (is.null(getGeneric("filterInjectionTime"))) {
    setGeneric("filterInjectionTime",
        function(object, ...)standardGeneric("filterInjectionTime"))
}
if (is.null(getGeneric("filterIntensity"))) {
    setGeneric("filterIntensity",
        function(object, ...)standardGeneric("filterIntensity"))
}
if (is.null(getGeneric("filterNonReplicatedFragments"))) {
    setGeneric("filterNonReplicatedFragments",
        function(object, ...) standardGeneric("filterNonReplicatedFragments"))
}
if (is.null(getGeneric("fragmentationMap"))) {
    setGeneric("fragmentationMap",
        function(object, ...)standardGeneric("fragmentationMap"))
}
if (is.null(getGeneric("removeEmptyConditions"))) {
    setGeneric("removeEmptyConditions",
        function(object, ...)standardGeneric("removeEmptyConditions"))
}
if (is.null(getGeneric("updateMedianInjectionTime"))) {
    setGeneric("updateMedianInjectionTime",
        function(object, ...)standardGeneric("updateMedianInjectionTime"))
}

