rTParam <- function() {
  # Constructor method for rTParam
  
  rTParam<-data.frame(
  ### This list all the parameters officially supported in the API
                      "list path, default parameters" = NA,       
                      "list path, taxonomy information" = NA,     
                      "output, histogram column width" = NA, 
                      "output, histograms" = NA,
                      "output, log path" = NA,
                      "output, maximum valid expectation value" = NA,
                      "output, message" = NA,
                      "output, one sequence copy" = NA,
                      "output, parameters" = NA,
                      "output, path" = NA,
                      "output, path hashing" = NA,
                      "output, performance" = NA,
                      "output, proteins" = NA,
                      "output, results" = NA,
                      "output, sequence path" = NA,
                      "output, sort results by" = NA,
                      "output, sequences" = NA,
                      "output, spectra" = NA,
                      "output, xsl path" = NA,
                      "protein, cleavage C-terminal mass change" = NA,
                      "protein, cleavage N-terminal mass change" = NA,
                      "protein, cleavage semi" = NA,
                      "protein, cleavage site" = NA,
                      "protein, C-terminal residue modification mass" = NA,
                      "protein, N-terminal residue modification mass" = NA,
                      "protein, modified residue mass file" = NA,
                      "protein, quick acetyl" = NA,
                      "protein, quick pyrolidone" = NA,
                      "protein, stP bias" = NA,
                      "protein, saps" = NA,
                      "protein, taxon" = NA,
                      "protein, use annotations" = NA,
                      "refine, cleavage semi" = NA,
                      "refine, maximum valid expectation value" = NA,
                      "refine, modification mass" = NA,
                      "refine, point mutations" = NA,
                      "refine, potential modification mass" = NA,
                      "refine, potential modification motif" = NA,
                      "refine, potential N-terminus modifications" = NA,
                      "refine, potential C-terminus modifications" = NA,
                      "refine, refine" = NA,
                      "refine" = NA,
                      "refine, saps" = NA,
                      "refine, sequence path" = NA,
                      "refine, spectrum synthesis" = NA,
                      "refine, tic percent" = NA,
                      "refine, unanticipated cleavage" = NA,
                      "refine, use annotations" = NA,
                      "refine, use potential modifications for full refinement" = NA,
                      "residue, modification mass" = NA,
                      "residue, potential modification mass" = NA,
                      "residue, potential modification motif" = NA,
                      "scoring, a ions" = NA,
                      "scoring, b ions" = NA,
                      "scoring, c ions" = NA,
                      "scoring, cyclic permutation" = NA,
                      "scoring, include reverse" = NA,
                      "scoring, maximum missed cleavage sites" = NA,
                      "scoring, minimum ion count" = NA,
                      "scoring, x ions" = NA,
                      "scoring, y ions" = NA,
                      "scoring, z ions" = NA,
                      "spectrum, contrast angle" = NA,
                      "spectrum, dynamic range" = NA,
                      "spectrum, fragment mass error" = NA,
                      "spectrum, fragment mass error units" = NA,
                      "spectrum, fragment mass type" = NA,
                      "spectrum, fragment monoisotopic mass error" = NA,
                      "spectrum, fragment monoisotopic mass error units" = NA,
                      "spectrum, minimum fragment mz" = NA,
                      "spectrum, minimum peaks" = NA,
                      "spectrum, minimum parent m+h" = NA,
                      "spectrum, neutral loss mass" = NA,
                      "spectrum, neutral loss window" = NA,
                      "spectrum, parent monoisotopic mass error minus" = NA,
                      "spectrum, parent monoisotopic mass error plus" = NA,
                      "spectrum, parent monoisotopic mass error units" = NA,
                      "spectrum, parent monoisotopic mass isotope error" = NA,
                      "spectrum, path" = NA,
                      "spectrum, path type" = NA,
                      "spectrum, sequence batch size" = NA,
                      "spectrum, skyline path" = NA,
                      "spectrum, threads" = NA,
                      "spectrum, total peaks" = NA,
                      "spectrum, use neutral loss window" = NA,
                      "spectrum, use noise suppression" = NA,
                      "spectrum, use contrast angle" = NA,

### This list the parameters that are not officially supported in the API, but are listed as implemented in the TPP documentation.                               
                      "output, http" = NA,
                      "output, sort best scores by" = NA,
                      "output, title" = NA,
                      "protein, cleavage N-terminal limit" = NA,
                      "protein, homolog management" = NA,
                      "protein, use minimal annotations" = NA,
                      "refine, full unanticipated cleavage" = NA,
                      "refine, potential N-terminus modification position limit" = NA,
                      "residue, NG deamidation" = NA,
                      "scoring, algorithm" = NA,
                      "scoring, pluggable scoring" = NA,
                      "spectrum, allowed neutral losses" = NA,
                      "spectrum, check all charges" = NA,
                      "spectrum, homolog error" = NA,
                      "spectrum, maximum parent charge" = NA,
                      "spectrum, use conditioning" = NA,
                      check.names=FALSE
                      )
  class(rTParam)<-c("data.frame", "rTParam")
  return(rTParam)
}

setParamValue <- function(param, category, parameter=NULL, value) {
  if( is.null(parameter) ){
    key <- category
  } else {
    key <- paste(category, parameter, sep=", ")
  }
  
  if (is.null(param[[key]])){
    warning("This category/parameter combination is not recognized by rTANDEM. It might be due to a typo, or to the use of an undocumented parameter. After the analysis, check your result's 'used.parameters' and 'unused.parameters' slots to confirm that the parameter was successfully used.")
  }
  param[[key]] <- value
  param
}
      
setParamDefault <-  function(param=NULL) {
  # Sets some default parameters values
  # Args:
  #    param: an optional rTParam object to be modified
  # Return:
  #    A rTParam object with parameters appropriate for an
  #    orbitrap instrument.

  param <- .checkParam(param)
  xls.path <- system.file("extdata/tandem-style.xsl",package='rTANDEM')
  myValues <- matrix(ncol=3, byrow=TRUE, data=list(
    "spectrum", "dynamic range", 100.0,
    "spectrum", "total peaks", 50,
    "spectrum", "use noise suppression", "no",
    "spectrum", "threads", 1,
    "spectrum", "minimum parent m+h", 500.0,
    "spectrum", "minimum fragment mz", 150.0,
    "spectrum", "minimum peaks", 15,
    "spectrum", "maximum parent charge", 4,
    "residue", "modification mass", "57.021464@C",
    "residue", "potential modification mass", "15.994915@M",
    "protein", "cleavage site", "[RK]|{P}",
    "protein", "N-terminal residue modification mass", 0.0,
    "protein", "C-terminal residue modification mass", 0.0,
    "scoring", "minimum ion count", 4,
    "scoring", "maximum missed cleavage sites", 1,
    "refine", NULL, "yes",
    "refine", "spectrum synthesis", "yes",
    "refine", "maximum valid expectation value", 0.1,
    "refine", "potential N-terminus modifications", NA,
    "refine", "potential C-terminus modifications", NA,
    "refine", "unanticipated cleavage", "yes",
    "refine", "potential modification mass", "15.994915@M,0.9848@N,0.9848@Q",
    "refine", "use potential modifications for full refinement", "no",
    "refine", "point mutations", "no",
    "output", "sort results by", "protein",
    "output", "path hashing", "no",
    "output", "xsl path", xls.path,
    "output", "parameters", "yes",
    "output", "performance", "yes",
    "output", "spectra", "yes",
    "output", "histograms", "yes",
    "output", "proteins", "yes",
    "output", "sequences", "yes",
    "output", "results", "valid",
    "output", "maximum valid expectation value", 0.1,
    "output", "histogram column width",	30) )

  for( i in 1:nrow(myValues) ){
    param <- setParamValue(param, myValues[[i,1]], myValues[[i,2]], myValues[[i,3]])
  }
  return(param)
}

setParamOrbitrap <- function(param=NULL) {
  # Sets some defaults parameters for orbitrap spectrometers.
  # Args:
  #    param: an optional rTParam object to be modified
  # Return:
  #    A rTParam object with parameters appropriate for an
  #    orbitrap instrument.
  param <- .checkParam(param)
  myValues <- matrix(ncol=3, byrow=TRUE, data=list(
    "spectrum", "fragment monoisotopic mass error", 0.4,
    "spectrum", "fragment monoisotopic mass error units", "Daltons",
    "spectrum", "parent monoisotopic mass error plus", 20,
    "spectrum", "parent monoisotopic mass error minus", 20,
    "spectrum", "parent monoisotopic mass error units", "ppm",
    "spectrum", "parent monoisotopic mass isotope error", "yes" ) )

  for( i in 1:nrow(myValues) ){
    param <- setParamValue(param, myValues[[i,1]], myValues[[i,2]], myValues[[i,3]])
  }
  return(param)
}

setParamQuadTof05da <- function(param=NULL) {
  # Sets some defaults parameters for XXXX spectrometers.
  # Args:
  #    param: an optional rTParam object to be modified
  # Return:
  #    A rTParam object with parameters appropriate for an
  #    XXXXXXXXX instrument.
  param <- .checkParam(param)
  myValues <- matrix(ncol=3, byrow=TRUE, data=list(
    "spectrum", "fragment monoisotopic mass error", 0.4,
    "spectrum", "fragment monoisotopic mass error units", "Daltons",
    "spectrum", "parent monoisotopic mass error plus", 0.5,
    "spectrum", "parent monoisotopic mass error minus",	0.5,
    "spectrum", "parent monoisotopic mass error units",	"Daltons",
    "spectrum", "parent monoisotopic mass isotope error", "no" ) )

  for( i in 1:nrow(myValues) ){
    param <- setParamValue(param, myValues[[i,1]], myValues[[i,2]], myValues[[i,3]])
  }
  return(param)
}

setParamQuadTof100ppm <- function(param=NULL) {
  # Sets some defaults parameters for XXXXXXXXX spectrometers.
  # Args:
  #    param: an optional rTParam object to be modified
  # Return:
  #    A rTParam object with parameters appropriate for an
  #    XXXXXXXXXX instrument.
  param <- .checkParam(param)

  myValues <- matrix(ncol=3, byrow=TRUE, data=list(
    "spectrum", "fragment monoisotopic mass error", 200,
    "spectrum", "fragment monoisotopic mass error units", "ppm",
    "spectrum", "parent monoisotopic mass error plus", 100,
    "spectrum", "parent monoisotopic mass error minus", 100,
    "spectrum", "parent monoisotopic mass error units",	"ppm",
    "spectrum", "parent monoisotopic mass isotope error", "yes" ) )

  for( i in 1:nrow(myValues) ){
    param <- setParamValue(param, myValues[[i,1]], myValues[[i,2]], myValues[[i,3]])
  }
  return(param)  
}

setParamIonTrap <- function(param=NULL) {
  # Sets some defaults parameters for ion trap spectrometers.
  # Args:
  #    param: an optional rTParam object to be modified
  # Return:
  #    A rTParam object with parameters appropriate for an
  #    ion trap instrument.
  param <- .checkParam(param)

  myValues <- matrix(ncol=3, byrow=TRUE, data=list(
    "spectrum", "fragment monoisotopic mass error", 0.4,
    "spectrum", "fragment monoisotopic mass error units", "Daltons",
    "spectrum", "parent monoisotopic mass error plus", 3.0,
    "spectrum", "parent monoisotopic mass error minus", 0.5,
    "spectrum", "parent monoisotopic mass error units",	"Daltons",
    "spectrum", "parent monoisotopic mass isotope error", "no" ) )
  
  for( i in 1:nrow(myValues) ){
    param <- setParamValue(param, myValues[[i,1]], myValues[[i,2]], myValues[[i,3]])
  }
  return(param)
}

rTTaxo <- function(taxon, format="peptide", URL) {
  # Constructor method for rTTaxo
  # Args:
  #    taxon : A taxon for the taxonomy (eg. "yeast") or a vector of taxa.
  #    format: The format of the database (eg. "peptides" or "spectrum") or a vector of formats.
  #    URL   : The path to the database file or a vector of paths.
  # Returns:
  #    A rTTaxo object.

  if( ! format %in% c("peptide", "spectrum", "saps", "mods") ){
    warning(paste("\"", format, "\" might not be recognized. The four formats of database for tandem are: peptide, spectrum, mods, and saps.", sep=""))
  }
  
  rTTaxo <- data.frame(
                       row.names=NULL,
                       taxon=taxon,
                       format=format,
                       URL=URL
                       )
  class(rTTaxo) <- c('rTTaxo', 'data.frame')
  return(rTTaxo)
}

addTaxon <- function(taxonomy=NULL, taxon, format="peptide", URL){
  # Adds a new taxon to a rTTaxo object.
  # Args:
  #    taxonomy: a rTTaxo object to which a new taxon will be added.
  #    taxon   : a string identifying the new taxon 
  #    format  : The format of the database ["peptide" | "saps" | "mods" | "spectrum" ]
  #    URLs    : a string or a vector containing the paths to the databases.
  # Returns:
  #    A rTTaxo object.
  if( is.null(taxonomy) ){
    return(rTTaxo(taxon=taxon, format=format,URL=URL))
  }
  if(! "rTTaxo" %in% class(taxonomy) ){
    stop("The taxonomy object must be of the class 'rTTaxo'.")
  }
  new.taxon <- rTTaxo(taxon=taxon, format=format, URL=URL)
  rTTaxo <- rbind(taxonomy, new.taxon)
  class(rTTaxo) <- c('rTTaxo', 'data.frame')
  return(rTTaxo)
}

.checkParam <- function(param) {
  # Used to check the validity of the 'param' object
  # passed to the various setParam functions. 
  if (! is.null(param) && ! "rTParam" %in% class(param) ) {
    stop("The parameter object must be of the class 'rTParam'")
  }
  if (is.null(param) ){ param <- rTParam() }
  return(param)
}

print.rTParam <- function(x, ...) {
  # Where x is a rTParam object.
  # print the defined parameter with minimal formating.
  cat(paste("rTParam object with", sum(!is.na(x)),"defined parameters.\n\n"))
  for(i in 1:length(x))
    if (! is.na(x[[i]]) ){
      cat(paste(names(x)[[i]],": \n\t", x[[i]], "\n\n", sep=""))
    }
}
