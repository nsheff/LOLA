# This is just a private test analysis for me to run the vignette tests
# pointing to local folders.
# Should probably put this elsewhere eventually.

library("LOLA")

dir = "~/linkto/resources/"

readRegionSetAnnotation(paste0(dir, "LOLACore/hg19"), refreshCaches=TRUE)
regionDB = loadRegionDB(paste0(dir, "LOLACore/hg19"))

regionSetA = readBed(paste0(dir, "LOLACore_deploy/lola_vignette_data/setA_100.bed"))
regionSetB = readBed(paste0(dir, "LOLACore_deploy/lola_vignette_data/setB_100.bed"))
regionSetC = readBed(paste0(dir, "LOLACore_deploy/lola_vignette_data/setC_100.bed"))
activeDHS = readBed(paste0(dir, "LOLACore_deploy/lola_vignette_data/activeDHS_universe.bed"))


userSets = GRangesList(regionSetA, regionSetB, regionSetC)
locResults = runLOLA(userSets, activeDHS, regionDB, cores=1)

userSets2 = GRangesList(regionSetB, regionSetC)
restUniv = buildRestrictedUniverse(userSets2)
locResults2 = runLOLA(userSets2, restUniv, regionDB, cores=1)

locResults[1:10,]
locResults2[1:10,]

locResults[order(maxRnk),][userSet==1,][1:10,]
locResults[order(maxRnk),][userSet==2,][1:10,]
locResults[order(maxRnk),][userSet==3,][1:10,]



locResults2[order(maxRnk),][userSet==1,][1:25,]
locResults2[order(maxRnk),][userSet==2,][1:25,]

writeCombinedEnrichment(locResults2, "out")


