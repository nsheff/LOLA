############################################################################################## 
# PEP-loading function 

##############################################################################################


#' Helper function that allows LOLA to read PEP-formatted databases 
#' 
#' @param config_location folder where the project_config.yaml file is located
#' @return PEPdB contains database location, samples annotation sheet, .yaml config file 
#' and the regions Granges list 
#' @export
#' @examples 
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' PEPdb = loadPEPdb(configLocation=dbpath)
#' 


loadPEPdb = function(configLocation){
  configLoc = file.path(path=configLocation)  
  if (file.exists(configLoc)) {
      pep_obj = pepr::Project(file = configLoc)
      # include samples annotation and config.yaml in list output
      samplesAnnotation = pepr::samples(pep_obj)
      configFile = pepr::config(pep_obj)
      samplesdf = as.data.frame(samplesAnnotation)
      chromRanges = lapply(samplesdf$file_path, LOLA::readBed) # need to make list into a GRanges list object
      regions = GRangesList(chromRanges)
      return(list(dblocation = configLoc,
                  configYAML = configFile,
                  regionAnno = samplesAnnotation[, -c("file_path", "genome")],
                  regionGRL = regions))
  } else {
      stop("could not find .yaml config file")
  }
}

#pep_object = loadPEPdb("/project/shefflab/resources/regions/LOLAHema/LOLAHema_PEP/project_config.yaml")




