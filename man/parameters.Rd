\name{parameters}
\alias{rTParam}
\alias{setParamValue}
\alias{SetParamValue}
\alias{setParamDefault}
\alias{SetParamDefault}
\alias{setParamPTMTreeSearch}
\alias{SetParamPTMTreeSearch}
\alias{setParamOrbitrap}
\alias{SetParamOrbitrap}
\alias{setParamIonTrap}
\alias{SetParamIonTrap}
\alias{setParamQuadTof05da}
\alias{SetParamQuadTof05da}
\alias{setParamQuadTof100ppm}
\alias{SetParamQuadTof100ppm}
\alias{rTTaxo}
\alias{addTaxon}
\alias{AddTaxon}
\alias{print.rTParam}
\title{
Creates and manipulates parameter objects for rTANDEM
}
\description{
Those functions instantiate parameter objects for rTANDEM and
manipulate their values. Many of those functions give default values,
either for general settings (e.g. setParamDefault), for instrument
specific settings (e.g. setParamOrbitrap), or for
function specific settings (e.g. setParamPTMTreeSearch).
}
\usage{
rTParam()
setParamValue(param, category, parameter=NULL, value)
setParamDefault(param=NULL)
setParamPTMTreeSearch(param=NULL)
setParamOrbitrap(param=NULL)
setParamIonTrap(param=NULL)
setParamQuadTof05da(param=NULL)
setParamQuadTof100ppm(param=NULL)
rTTaxo(taxon=NA, format=NA, URL=NA)
addTaxon(taxonomy=NULL, taxon, format="peptide", URL)
\method{print}{rTParam}(x, ...)
}
\arguments{
  \item{param}{
    A parameter object of class rTParam.
  }
  \item{x}{
    A parameter object of class rTParam.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
  \item{category}{
    A string representing a category of parameters (e.g. "output",
    "protein").
  }
  \item{parameter}{
    A string representing the name of a parameter.
  }
  \item{value}{
    The value to give to the parameter.
  }
  \item{taxon}{
    A string representing the name of a taxon (e.g. "homo sapiens").
  }
  \item{taxonomy}{
    A taxonomy object of class rTTaxo.
  }
  \item{format}{
    A string representing the type of the database. The four types are:
    "peptide", "saps", "mods" and "spectrum".
  }
  \item{URL}{
    A string representing the full path to the database file.
  }
}
\value{
  The functions rTParam, and all functions setParam..., return an
  object of S3 class rTParam.
  The functions rTTaxo and addTaxon return an object of S3 class rTTaxo.
  print.rTParam is used for its side effect of displaying a parameter object.
}
\references{
  \url{thegpm.org/TANDEM/api/}
}
\seealso{
  \code{\link{GetParamFromXML}}, \code{\link{GetTaxoFromXML}},
  \code{\link{WriteTaxoToXML}}, \code{\link{WriteParamToXML}}.
}
\examples{
# Initialize an empty parameter object:
param <- rTParam()
print.rTParam(param)

# Set general values.
param <- setParamDefault(param)

# Add instrument specific values for a LTQ mass spectrometer.
param <- setParamIonTrap(param)
print.rTParam(param)
}
