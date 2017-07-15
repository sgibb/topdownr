setClass("TopDownExperiment",
  slots=c(name="character",
          sequence="character",
          fragmentTable="data.table",
          assignmentTable="data.table"),
  contains=c("MSnExp"),
  prototype=prototype(new("MSnExp"),
                      name=character(),
                      sequence=character(),
                      fragmentTable=data.table(),
                      assignmentTable=data.table()),
  validity=function(object).validateTopDownExperiment(object)
)

## just subclass to ensure that the elementMetadata slot (DataFrame) contains
## a mass column
setClass("FragmentViews",
         contains="XStringViews",
         validity=function(object).validateFragmentViews(object))

setClass("TopDownSet",
  slots=c(rowViews="FragmentViews",
          colData="DataFrame",
          assays="CsparseMatrix"),
  prototype=prototype(rowViews=new("FragmentViews"),
                      colData=new("DataFrame"),
                      assays=new("dgCMatrix"))
  #validity=function(object).validateTopDownSet(object)
)
