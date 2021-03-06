setClass("rTResult",
         representation=representation(
           result.file                        ="character",
           proteins                           ="data.table",
           peptides                           ="data.table",
           ptm                                ="data.table",
           used.parameters                    ="data.frame",
           unused.parameters                  ="vector",
           sequence.source.paths              ="vector",
           estimated.false.positive           ="integer",
           total.peptides.used                ="integer",
           total.proteins.used                ="integer",
           total.spectra.assigned             ="integer",
           total.spectra.used                 ="integer",
           total.unique.assigned              ="integer",
           start.time                         ="character",
           xtandem.version                    ="character",
           quality.values                     ="vector",
           nb.input.models                    ="integer",
           nb.input.spectra                   ="integer",
           nb.partial.cleavages               ="integer",
           nb.point.mutations                 ="integer",
           nb.potential.C.terminii            ="integer",
           nb.potential.N.terminii            ="integer",
           nb.unanticipated.cleavages         ="integer",
           initial.modelling.time.total       ="numeric",
           initial.modelling.time.per.spectrum="numeric",
           load.sequence.models.time          ="numeric",
           refinement.per.spectrum.time       ="numeric"
           )
         )

# rTResult_s is an extension of rTResult that accomodates spectra.
setClass("rTResult_s",
         contains="rTResult",
         representation=representation(spectra="data.table")
         )

setIs("rTResult_s", "rTResult")
