## just subclass to ensure that the elementMetadata slot (DataFrame) contains
## a mass column
setClass("FragmentViews",
         contains="XStringViews",
         validity=function(object).validateFragmentViews(object))

setClass("TopDownSet",
  slots=c(rowViews="FragmentViews",
          colData="DataFrame",
          assay="Matrix",
          files="character",
          tolerance="numeric",
          processing="character"),
  prototype=prototype(rowViews=new("FragmentViews"),
                      colData=new("DataFrame"),
                      assay=new("dgCMatrix"),
                      files=character(),
                      tolerance=double(),
                      processing=character()),
  validity=function(object).validateTopDownSet(object)
)
