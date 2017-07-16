## just subclass to ensure that the elementMetadata slot (DataFrame) contains
## a mass column
setClass("FragmentViews",
         contains="XStringViews",
         validity=function(object).validateFragmentViews(object))

setClass("TopDownSet",
  slots=c(rowViews="FragmentViews",
          colData="DataFrame",
          assays="Matrix",
          files="character",
          processing="character"),
  prototype=prototype(rowViews=new("FragmentViews"),
                      colData=new("DataFrame"),
                      assays=new("dgCMatrix"),
                      files=character(),
                      processing=character()),
  validity=function(object).validateTopDownSet(object)
)
