############################################################################################## 
# PEP-loading function 

##############################################################################################


#' Helper function that allows LOLA to read PEP-formatted databases 
#' 
#' @param configFolder     folder where the config.yaml file is located
#' @param configName       Name of the config.yaml file 
#' @param useCache            uses simpleCache to cache and load the results    
#' @return PEPdB contains database location, samples annotation sheet, .yaml config file 
#' and the regions Granges list 
#' @export
#' @examples 
#' ConfigPath = system.file("extdata", "hg19", "ucsc_examplePEP", "project_config.yaml", package="LOLA")
#' PEPdb = loadPEPdb_CACHE(configLocation=ConfigPath)
#' 


loadPEPdb_CACHE = function(configFolder, configName, useCache=TRUE){
    configLoc = file.path(path=paste0(configFolder, configName))
    if (file.exists(configLoc)) {
        # Use pepr to read in the PEP metadata
        pepObject = pepr::Project(file = configLoc)
        configFile = pepr::config(pepObject)
        samplesAnnotation = pepr::samples(pepObject)
        # Output the samples annotation and make it a dataframe so that we can lapply through each file path
        samplesdf = as.data.frame(samplesAnnotation)
        if (useCache & requireNamespace("simpleCache", quietly=TRUE)){
            simpleCache::simpleCache("chromRanges", { # need to make regions GRanges objects and cache the data
                lapply(samplesdf$file_path, LOLA::readBed)}, 
                                     cacheDir=file.path(path=configFolder),
                                     recreate=FALSE) 
         } else {
             if (nrow(samplesAnnotation) > 100) {
                 # tell the user they should install simpleCache 
                 message("You should istall simpleCache so that you save time when loading your database next")
          }
             chromRanges = lapply(samplesdf$file_path, LOLA::readBed)    
      }
      regions = GRangesList(chromRanges)
      return(list(configLocation = configLoc,
                  configYAML = configFile,
                  regionAnno = samplesAnnotation[, -c("genome")],
                  regionGRL = regions))
    } else {
        stop("could not find .yaml config file")
  }
}









