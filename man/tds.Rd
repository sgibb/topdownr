% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tds.R
\docType{data}
\name{tds}
\alias{tds}
\title{TopDownSet Example Data}
\format{
A \linkS4class{TopDownSet}
with 14901 fragments (1949 rows, 351 columns).
}
\source{
Subset taken from the
\link[topdownrdata:topdownrdata-package]{topdownrdata::topdownrdata-package}
package.
}
\usage{
tds
}
\description{
An example data set for \code{topdownr}. It is just a subset of the myoglobin
dataset available in
\link[topdownrdata:topdownrdata-package]{topdownrdata::topdownrdata-package}.
}
\details{
It was created as
follows:

\if{html}{\out{<div class="sourceCode">}}\preformatted{tds <- readTopDownFiles(
   topdownrdata::topDownDataPath("myoglobin"),
   ## Use an artifical pattern to load just the fasta
   ## file and files from m/z == 1211, ETD reagent
   ## target 1e6 and first replicate to keep runtime
   ## of the example short
   pattern=".*fasta.gz$|1211_.*1e6_1",
   adducts=data.frame(mass=1.008, name="zpH", to="z"),
   neutralLoss=PSMatch::defaultNeutralLoss(
       disableWaterLoss=c("Cterm", "D", "E", "S", "T")),
   tolerance=25e-6)
}\if{html}{\out{</div>}}
}
\keyword{datasets}
