# Unit tests
library(LOLA)

context("Testthat context...")

test_that("loadRegionDB",  {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	expect_equal(length(regionDB$regionGRL), 5)
	expect_equal(length(unlist(regionDB$regionGRL)), 32749)
	expect_identical(regionDB$regionGRL[[4]], regionDB$regionGRL[[5]])
	expect_identical(regionDB$collectionAnno$collector, "Nathan")
	expect_identical(regionDB$collectionAnno$source, "UCSC Genome Browser")

})

test_that( "runLOLA", {
	dbPath = system.file("extdata", "hg19", package="LOLA")
	regionDB = loadRegionDB(dbPath)
	data("sample_input", package="LOLA") # load userSet
	data("sample_universe", package="LOLA") # load userUniverse
	locResults = runLOLA(userSet, userUniverse, regionDB, cores=1)

	expect_equal(nrow(locResults), 5)
	expect_true(all(locResults[,support] == c(662, 121, 121, 4, 3006)))
	expect_true(all(locResults[,filename] ==
	c("laminB1Lads.bed", "vistaEnhancers.bed", "vistaEnhancers_colNames.bed", "numtSAssembled.bed", "cpgIslandExt.bed")))
	expect_equal(nrow(locResults), 5)


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
# To do: write a test case for splitFileIntoCollection():
#splitFileIntoCollection(system.file("extdata", "examples/combined_regions.bed", package="LOLA"))
test_that("readBed", {
	file = system.file("extdata", "examples/combined_regions.bed", package="LOLA")
	rb = readBed(file)
	expect_equal(length(rb), 16)
	# Make sure strand is getting picked up correctly:
	expect_false("*" %in% as.character(strand(rb)))
})
