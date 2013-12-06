ms2.plot <- function(spectrum.id, result) {
  ## Plot a raw ms2 spectra
  ## Args:
  ##   spectrum.id: The id of a spectrum (from the field result@spectra$id)
  ##   result: A result object of class rTResult_s
  ## Returns:
  ##   A plot of the ms2 spectra

  ## To prevent "no visible binding" warning with data.table accessors
  id=spectrum.mh=spectrum.maxI=NULL
  rm(id, spectrum.mh, spectrum.maxI)
  
  if( ! is(result, "rTResult_s") ){
    stop("The 'result' parameter is of class: \"",
         class(result),
         "\" it must be of class 'rTResult_s'.")
  }

  if ( ! spectrum.id %in% result@spectra$id ){
    stop("There is no spectrum \"", spectrum.id,
         "\" in the given dataset.")
  }
  
  spectra  <- as.list(subset(result@spectra, id==spectrum.id))
  pep.data <- subset(result@peptides,
                     spectrum.id==spectrum.id,
                     select=c(spectrum.mh, spectrum.maxI))[1]
  mh      <- pep.data[[1]]
  maxI    <- pep.data[[2]]
  Xdata   <- as.numeric(strsplit(spectra$Xdata, "\\s+")[[1]])
  Ydata   <- as.numeric(strsplit(spectra$Ydata, "\\s+")[[1]])
  Xunit   <- spectra$Xunit
  Yunit   <- spectra$Yunit
  if(Xunit=="MASSTOCHARGERATIO"){Xunit="Mass To Charge Ratio"}
  if(Yunit=="UNKNOWN") {Yunit="Normalized intensity"}
  plot(Xdata,Ydata, type="h",
       xlab=Xunit, ylab=Yunit,
       main=paste("Spectrum: \"",spectra$label,"\"\n\nMaximum intensity: ", maxI, "     mh: ", mh, sep=""))
}
