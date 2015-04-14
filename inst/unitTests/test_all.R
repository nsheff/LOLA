# Unit tests

test_readRegionSetAnnotation = function() {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionAnno = readRegionSetAnnotation(dbLocation= dbPath)
	checkEquals(	nrow(regionAnno), 4)
}

test_loadRegionDB = function() {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	checkEquals(length(regionDB$regionGRL), 4)
	checkEquals(length(unlist(regionDB$regionGRL)), 31410)
}

test_calcLocEnrichment = function() {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	data("sample_input", package="LOLA") # load userSet
	data("sample_universe", package="LOLA") # load userUniverse
	locResults = calcLocEnrichment(userSet, userUniverse, regionDB, dbTitle="dbTitle", cores=1)

	checkEquals(nrow(locResults), 4)
	checkTrue(all(locResults[,support] == c(662, 121, 4, 3006)))
	checkTrue(all(locResults[,filename] == 
	c("laminB1Lads.bed", "vistaEnhancers.bed", "numtSAssembled.bed", "cpgIslandExt.bed")))
}

