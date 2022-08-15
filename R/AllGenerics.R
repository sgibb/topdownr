setGeneric("assayData",
    function(object)standardGeneric("assayData"))
setGeneric("rowViews",
    function(object, ...)standardGeneric("rowViews"))
setGeneric("colData",
    function(object, ...)standardGeneric("colData"))
setGeneric("colData<-",
    function(object, ..., value)standardGeneric("colData<-"))
setGeneric("conditionData",
    function(object, ...)standardGeneric("conditionData"))
setGeneric("conditionData<-",
    function(object, ..., value)standardGeneric("conditionData<-"))
setGeneric("conditionNames",
    function(object)standardGeneric("conditionNames"))
setGeneric("bestConditions",
    function(object, ...)standardGeneric("bestConditions"))
setGeneric("filterCv",
    function(object, ...)standardGeneric("filterCv"))
setGeneric("filterInjectionTime",
    function(object, ...)standardGeneric("filterInjectionTime"))
setGeneric("filterIntensity",
    function(object, ...)standardGeneric("filterIntensity"))
setGeneric("filterNonReplicatedFragments",
    function(object, ...) standardGeneric("filterNonReplicatedFragments"))
setGeneric("fragmentationMap",
    function(object, ...)standardGeneric("fragmentationMap"))
setGeneric("removeEmptyConditions",
    function(object, ...)standardGeneric("removeEmptyConditions"))
setGeneric("updateConditionNames",
    function(object, ...)standardGeneric("updateConditionNames"))
setGeneric("updateMedianInjectionTime",
    function(object, ...)standardGeneric("updateMedianInjectionTime"))
