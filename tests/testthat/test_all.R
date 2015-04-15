# Unit tests
library(LOLA)

context("Context here...")

test_that("loadRegionDB",  {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	expect_equal(length(regionDB$regionGRL), 5)
	expect_equal(length(unlist(regionDB$regionGRL)), 32749)
	expect_identical(regionDB$regionGRL[[4]], regionDB$regionGRL[[5]])
})

test_that( "calcLocEnrichment", {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	data("sample_input", package="LOLA") # load userSet
	data("sample_universe", package="LOLA") # load userUniverse
	locResults = calcLocEnrichment(userSet, userUniverse, regionDB, dbTitle="dbTitle", cores=1)

	expect_equal(nrow(locResults), 5)
	expect_true(all(locResults[,support] == c(662, 121, 121, 4, 3006)))
	expect_true(all(locResults[,filename] == 
	c("laminB1Lads.bed", "vistaEnhancers.bed", "vistaEnhancers_colNames.bed", "numtSAssembled.bed", "cpgIslandExt.bed")))
})

test_that("readRegionSetAnnotation", {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionAnno = readRegionSetAnnotation(dbLocation= dbPath)
	expect_equal(	nrow(regionAnno), 5)
})


