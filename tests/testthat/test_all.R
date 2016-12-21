# Unit tests
library(LOLA)

context("Testthat context...")

test_that("loadRegionDB",  {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	regionDB_col = loadRegionDB(dbPath, collections="ucsc_example")
		regionDB_col
	expect_identical(regionDB, regionDB_col)

	expect_error(loadRegionDB(paste(dbPath, "badfolder"), collections="ucsc_example"))


	expect_equal(length(regionDB$regionGRL), 5)
	expect_equal(length(unlist(regionDB$regionGRL)), 32749)
	expect_identical(regionDB$regionGRL[[4]], regionDB$regionGRL[[5]])
	expect_identical(regionDB$collectionAnno$collector, "Nathan")
	expect_identical(regionDB$collectionAnno$source, "UCSC Genome Browser")

	data("sample_input", package="LOLA") # load userSet
	data("sample_universe", package="LOLA") # load userUniverse
	userSet = userSets[[1]]
	# Test redefined user sets:
	userSetsRedefined =	redefineUserSets(list(userSet), userUniverse)
	fo = findOverlaps(userSetsRedefined[[1]], userUniverse, type='equal')
	expect_equal(length(fo), 3043)

	# Test restricted universe:
	ru = buildRestrictedUniverse(GRangesList(userSet[1:50], userSet[51:100]))
	expect_equal(length(ru), 100)

	dbPathMulti = system.file("extdata", "multi", package="LOLA")
	dbPathEmpty = system.file("extdata", "empty", package="LOLA")

	expect_equal(length(listRegionSets(dbPath)), 5)
	expect_equal(length(listRegionSets(dbPathMulti)), 6)

	expect_equal(sapply(getRegionSet(dbPath, c("cpgIslandExt.bed", "vistaEnhancers.bed")), length), c(28691, 1339))

	expect_error(loadRegionDB(dbPathEmpty))

	# check auto annotation file:
	rdb = loadRegionDB(dbPathMulti, collections=c("collection2"))
	expect_equal(rdb$collectionAnno$collectionname, "collection2")
	mergedDB = loadRegionDB(c(dbPathMulti, dbPath))
	expect_equal(nrow(mergedDB$collectionAnno), 3)
})

test_that( "runLOLA", {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	data("sample_input", package="LOLA") # load userSet
	userSet = userSets[[1]]
	data("sample_universe", package="LOLA") # load userUniverse
	locResults = runLOLA(userSet, userUniverse, regionDB, cores=1)


	expect_equal(nrow(locResults), 5)
	expect_true(all(locResults[,support] == c(632, 124, 124, 8, 3002)))
	expect_true(all(locResults[,filename] ==
	c("laminB1Lads.bed", "vistaEnhancers.bed", "vistaEnhancers_colNames.bed", "numtSAssembled.bed", "cpgIslandExt.bed")))
	expect_equal(nrow(locResults), 5)

	# Test minoverlaps:
	locResults1500 = runLOLA(userSet, userUniverse, regionDB, cores=1, minOverlap=1500)
locResults1500[,support]	
	expect_true(all(locResults1500[,support] == c(358, 425, 73, 73, 8)))
	locResult = locResults[2,]
	# Test post-enrichment functions:
	eeo = extractEnrichmentOverlaps(locResult, userSet, regionDB)
	expect_equal(length(eeo), 307)

	# Test writing results:
	# Make sure checkUniverseAppropriateness can run:
	checkUniverseAppropriateness(userSets, userUniverse)
	checkUniverseAppropriateness(userSets, userUniverse, fast=TRUE)
	locResultsMult = runLOLA(userSets, userUniverse, regionDB, cores=1, minOverlap=1500)

	# check bad universe:
	expect_warning(runLOLA(userSets, userUniverse[1:100], regionDB, cores=1, minOverlap=1500), "Negative b")

	extData = system.file("extdata", package="LOLA")
	tmpFolder= paste0(extData, "/test_temp")
	writeCombinedEnrichment(locResultsMult, tmpFolder)

	locResultRead = fread(paste0(tmpFolder, "/userSet_setA.tsv"))
	expect_equal(nrow(locResultRead), 5)
	locResultRead = fread(paste0(tmpFolder, "/userSet_setB.tsv"))
	expect_equal(nrow(locResultRead), 5)

	unlink(tmpFolder, recursive=TRUE)
})

test_that("readRegionSetAnnotation", {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionAnno = readRegionSetAnnotation(dbLocation= dbPath)
	expect_equal(	nrow(regionAnno), 5)
})

test_that("GRangesOverlaps", {
	r1 = GRanges("chrX", IRanges(start=c(1, 5, 9, 13, 27, 30),
		end=c(3, 7, 11, 21, 29, 32)))
	r2 = GRanges("chrX", IRanges(start=c(4, 6, 17, 20, 24, 42),
		end=c(5, 15, 18, 22, 28, 45)))
	#r2=reduce(r2);r2

	#r2 = disjoin(r2)
	isDisjoint(r2)
	r1
	r2
	# Each r1 region, Number of r2 regions overlapping each r1 region.
	countOverlaps(r1, r2)

	# Each r1 region, Overlaps anything in r2? Yes/no.
	#countOverlapsAny(r1, r2)
	expect_equal(countOverlaps(r1, GRangesList(r2)), c(0,1,1,1,1,0)) # identical

	#Total # of regions in r2 overlapping anything in r1 (sum(countOverlapsAny(r2,r1))
	expect_equal(countOverlaps(GRangesList(r1, r1), r2), c(5,5)) #5

	#Total # of regions in r1 overlapping anything in r2 (sum(countOverlapsAny(r1,r2))
	expect_equal(countOverlaps(GRangesList(r2), r1), 4) #4


	# Each r2 region, Overlaps anything in r1? Yes/no.
	#countOverlapsAny(r2, r1)
	countOverlaps(r2, GRangesList(r1)) #identical; but faster.

	# Total Number of regions in each r1 SET overlapping anything in r2;
	#lapplyAlias( list(r1, r1), countOverlapsRev, GRangesList(r2));

})


context("Test reading functions")
test_that("readBed", {
	cr = system.file("extdata", "examples/combined_regions.bed", package="LOLA")
	rb = readBed(cr)
	expect_equal(length(rb), 16)
	# Make sure strand is getting picked up correctly:
	expect_false("*" %in% as.character(strand(rb)))
	s = splitFileIntoCollection(cr, 4)
	i = fread(paste0(cr, "_collection/insulator.bed"))
	p = fread(paste0(cr, "_collection/promoter.bed"))
	e = fread(paste0(cr, "_collection/enhancer.bed"))
	expect_equal(nrow(e), 5)
	expect_equal(nrow(p), 8)
	expect_equal(nrow(i), 3)

	# Make sure the bed reading is switching to 1-based coordinates
	expect_equal(start(ranges(rb))[1], 28736)
	unlink(paste0(cr, "_collection"), recursive=TRUE)
})

