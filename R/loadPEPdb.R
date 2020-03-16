############################################################################
# PEP-loading function 

############################################################################


#' Given a PEP-formatted database, this function will display its metadata and 
#' load the genomic regions for runLOLA() overlap calculations
#' @param configPath Psth to the PEP config.yaml file  
#' @param useCache   uses simpleCache to cache and load the results    
#' @return A list containing the config location, config file,   
#' samples annotation sheet and the regions Granges list 
#' @export
#' @examples 
#' ConfigPath = system.file("extdata/hg19/ucsc_examplePEP/project_config.yaml",
#'                         package="LOLA")
#' PEPdb = loadPEPdb(ConfigPath)
#' 



loadPEPdb = function(configPath, useCache=TRUE){
	configLoc = file.path(path=configPath)
	if (!file.exists(configLoc)) {
		stop("Could not find .yaml config file")
	} else {
		# Use pepr to read in the PEP metadata
		pepObject = pepr::Project(file = configLoc)
		configFile = pepr::config(pepObject)
		smpl = pepr::samples(pepObject)
		# Rename columns to make annotation LOLA compatible 
		colnames(smpl)[colnames(smpl)=="file_name"] = "filename"
		colnames(smpl)[colnames(smpl)=="cell_type"] = "cellType"
		colnames(smpl)[colnames(smpl)=="exp_protocol"] = "collection"
		colnames(smpl)[colnames(smpl)=="data_source"] = "dataSource"
		# Make samples annotation a df and lapply through each file path
		samplesdf = as.data.frame(smpl)
		configFolder = dirname(configPath)
		if (useCache & requireNamespace("simpleCache", quietly=TRUE)){
			simpleCache::simpleCache("chromRanges", { 
			lapply(samplesdf$output_file_path, LOLA::readBed)}, 
			cacheDir=file.path(path=configFolder), 
			recreate=FALSE) 
		} else {
			if (nrow(samplesAnnotation) > 100) {
				# tell the user they should install simpleCache 
				message("You should install simpleCache 
					to save time when you load the 
					database next")
			}
			chromRanges = lapply(samplesdf$output_file_path, LOLA::readBed)    
		}
		regions = GRangesList(chromRanges)
		igdDBlocation = configFile$iGD_dir
		return(list(configLocation = configLoc,
				configYAML = configFile,
				regionAnno = smpl[, -c("genome")],
				regionGRL = regions,
				iGDRefDatabase = igdDBlocation))
	}
}









