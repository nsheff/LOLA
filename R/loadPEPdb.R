############################################################################################## 
# PEP-loading function 

##############################################################################################


#' Helper function that allows LOLA to read PEP-formatted databases 
#' @param config_location folder where the project_config.yaml file is located 


loadPEPdb = function(...){
  config_loc = file.path(...) # feed config location into Project function 
  if (file.exists(config_loc)) {
    pep_obj = pepr::Project(file = config_loc)
    # include samples annotation and config.yaml in list output
    samples_annotation = samples(pep_obj)
    config_file = config(pep_obj)
    samples_df = as.data.frame(samples_annotation)
    chrom_ranges = lapply(samples_df$file_path, readBed) # need to make list into a granges list object
    regions = GRangesList(chrom_ranges)
    return(list(config_loc, samples_annotation,config_file, regions))
    
  } else {
    stop("could not find .yaml config file :(")
    }
}

#pep_object = readPEP("/project","shefflab","resources","regions","LOLAHema","PEP_project","project_config.yaml")


