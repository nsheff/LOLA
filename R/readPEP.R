############################################################################################## 
# PEP-loading function 

##############################################################################################


# Write a function that allows LOLA to use pepr so that it can read PEP-formatted databases 


readPEP = function(...){
  library(pepr)
  config = file.path(...)
  pep_database = pepr::Project(file = config)
  return(pep_database)
}

#pep_object = readPEP("/project","shefflab","resources","regions","LOLAHema","PEP_project","project_config.yaml")

#samples(pep_object)
#config(pep_object)
