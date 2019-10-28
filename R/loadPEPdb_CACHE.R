############################################################################################## 
# PEP-loading function 

##############################################################################################


#' Helper function that allows LOLA to read PEP-formatted databases 
#' 
#' @param config_location     folder where the project_config.yaml file is located
#' @param useCache            uses simpleCache to cache and load the results    
#' @return PEPdB contains database location, samples annotation sheet, .yaml config file 
#' and the regions Granges list 
#' @export
#' @examples 
#' ConfigPath = system.file("extdata", "hg19", "ucsc_examplePEP", "project_config.yaml", package="LOLA")
#' PEPdb = loadPEPdb_CACHE(configLocation=ConfigPath)
#' 


loadPEPdb_CACHE = function(configLocation, useCache=TRUE){
    configLoc = file.path(path=configLocation)  
    if (file.exists(configLoc)) {
        # Use pepr to read in the PEP metadata
        pep_obj = pepr::Project(file = configLoc)
        configFile = pepr::config(pep_obj)
        #Output the samples annotation and make it a dataframe so that we can lapply through each file path
        if (useCache & requireNamespace("simpleCache", quietly=TRUE)){
            setCacheDir("RCACHE.DIR")
            simpleCache::simpleCache("samplesAnnotation", {pepr::samples(pep_obj)}, 
                                     #cacheDir=RCACHE.dir,
                                     recreate=FALSE)  
            samplesdf = as.data.frame(samplesAnnotation)  
            simpleCache::simpleCache("chromRanges", {
                lapply(samplesdf$file_path, LOLA::readBed)},
                                     #cacheDir=RCACHE.dir,
                                     recreate=FALSE) # need to make list into a GRanges and cache the region outputs
         } else {
             samplesAnnotation = pepr::samples(pep_obj)
             chromRanges = lapply(samplesdf$file_path, LOLA::readBed)
             if (nrow(samplesAnnotation) > 100) {
                 # tell the user they should install simplecache 
                 message("You should istall simpleCache so that you save time next time you load your database")
        }
      }
      regions = GRangesList(chromRanges)
      return(list(dblocation = configLoc,
                configYAML = configFile,
                regionAnno = samplesAnnotation[, -c("genome")],
                regionGRL = regions))
    } else {
        stop("could not find .yaml config file")
  }
}









