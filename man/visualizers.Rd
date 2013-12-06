\name{visualizers}
\alias{ms2.plot}
\title{
Plot spectrum from a rTResult object
}
\description{
Functions like \code{GetTaxoFromXML("pathToXML")},
\code{GetParamFromXML("pathToXML")}, \code{GetResultsFromXML(pathToXML)}
creates R object from X!Tandem-style xml files. The functions
\code{WriteParamToXML(paramObject)} and
\code{WriteTaxoToXML(paramObject)} will create an X!Tandem-style xml
file from an R object.
}
\usage{
ms2.plot(spectrum.id, result)
}
\arguments{
  \item{spectrum.id}{
    The id of the spectrum to be plotted (from the field result@spectra$id).
  }
  \item{result}{
    A result object of class rTResult_s. 
  }
}
\value{
  'plot.ms2' returns a plot of the spectrum. 
}
\examples{

require(rTANDEM)
result <- GetResultsFromXML(
  system.file("extdata/result.xml", package="rTANDEM")
)

## Get the first spectra of the dataset and plot it:
spectrum.id <- result@spectra$id[[1]]
ms2.plot(spectrum.id, result)
}
