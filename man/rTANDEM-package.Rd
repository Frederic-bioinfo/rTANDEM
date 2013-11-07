\name{rTANDEM-package}
\alias{rTANDEM-package}
\docType{package}
\title{
An R encapsulation of X!TANDEM ('Jackhammer' release, 2013.06.15)
}
\description{
X!Tandem is an open source software for protein identification by
tandem mass spectrometry, and rTANDEM encapsultates this software in R.
rTANDEM provides a basic encapsulation of X!Tandem: it has a function
that takes as an argument the path to an X!Tandem style parameter file
and return the path to an X!tandem style output file. The package
also presents functions to transform parameters or results files
into R objects and vice versa. Some functions are also available to
examine the results within R. An associated package, shinyTANDEM,
provides a graphical interface to visualize results objects. 
}
\details{
\tabular{ll}{
Package: \tab rTANDEM\cr
Type: \tab Package\cr
Version: \tab 1.2.3\cr
Date: \tab 2013-11-07\cr
License: \tab Artistic License 1.0\cr
}
}
\author{
Authors: Frederic Fournier, Charles Joly Beauparlant, Rene Paradis,
Arnaud Droit
Maintainer: Frederic Fournier <frederic.fournier@crchuq.ulaval.ca>,
Charles Joly Beauparlant <charles.joly-beauparlant@crchul.ulaval.ca>
}
\references{
Robertson Craig and Ronald C. Beavis, TANDEM: matching proteins with
mass spectra, Bioinformatics, 2004, 20 1466-7.
http://www.thegpm.org/tandem/ 
}
\keyword{ package }
\examples{
# X!Tandem call style: we call tandem(input) on a single
# rTParam object.

# We create rTParam and from X!Tandem xml files
# located in the installation folder: 
param <- GetParamFromXML(system.file("extdata/input.xml", package="rTANDEM"))

# We create a rTTaxo object and identify a database for yeast
taxonomy <- rTTaxo( 
  taxon="yeast",
  format="peptide",
  URL= system.file("extdata/fasta/scd.fasta.pro", package="rTANDEM")
  )

# We will adjust those two objects to use one another and to use,
# the path of some data and configuration files located
# in the installation folder:
param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
param <- setParamValue(param, 'list path', 'default parameters',
  value=system.file("extdata/default_input.xml", package="rTANDEM"))
param <- setParamValue(param, 'spectrum', 'path',
  value=system.file("extdata/test_spectra.mgf", package="rTANDEM"))
param <- setParamValue(param, 'output', 'xsl path',
  value=system.file("extdata/tandem-input-style.xsl", package="rTANDEM"))
param <- setParamValue(param, 'output', 'path',
  value=paste(getwd(), "output.xml", sep="/"))

# This is the main command to run rTANDEM. The output will be
# written to a file in the working directory and the function
# returns the path to this file.
output.file <- tandem(param)
output.file
}
