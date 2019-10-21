############################################################################################## 
# PEP-loading function 

##############################################################################################


# Helper function that allows LOLA to read PEP-formatted databases 


readPEP = function(...){
  library(pepr)
  config_loc = file.path(...)
  if (file.exists(config_loc) & length(config_loc == 1)) {
    pep_obj = pepr::Project(file = config_loc)
    # include samples annotation and config.yaml in list output
    samples_annotation = samples(pep_obj)
    pep_config = config(pep_obj)
    samples_df = as.data.frame(samples_annotation)
    chrom_ranges = lapply(samples_df$file_path, fread) # need to make list into a granges list object
    regions = GRangesList(chrom_ranges)
    return(nlist(config_loc, samples_annotation, pep_config, regions))
    
  } else {
    stop("could not find .yaml config file")
    }
}

#pep_object = readPEP("/project","shefflab","resources","regions","LOLAHema","PEP_project","project_config.yaml")


