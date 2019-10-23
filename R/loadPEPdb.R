############################################################################################## 
# PEP-loading function 

##############################################################################################


#' Helper function that allows LOLA to read PEP-formatted databases 
#' 
#' @param config_location folder where the project_config.yaml file is located
#' @return PEPdB contains database location, samples annotation sheet, .yaml config file 
#' and the regions Granges list 


loadPEPdb = function(config_location){
  library(pepr)
  config_loc = file.path(path=config_location)  
  if (file.exists(config_loc)) {
      pep_obj = pepr::Project(file = config_loc)
      # include samples annotation and config.yaml in list output
      samples_annotation = samples(pep_obj)
      config_file = config(pep_obj)
      samples_df = as.data.frame(samples_annotation)
      chrom_ranges = lapply(samples_df$file_path, readBed) # need to make list into a GRanges list object
      regions = GRangesList(chrom_ranges)
      return(list(dblocation = config_loc, regionAnno = samples_annotation, regionGRL = regions))
  } else {
      stop("could not find .yaml config file")
  }
}

pep_object = loadPEPdb("/project/shefflab/resources/regions/LOLAHema/LOLAHema_PEP/project_config.yaml")


