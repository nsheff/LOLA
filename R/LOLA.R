# PACKAGE DOCUMENTATION
#' Genome locus overlap analysis.
#'
#' Run, Lola!
#'
#' @docType package
#' @name LOLA
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/sheffien}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#' @import BiocGenerics S4Vectors IRanges reshape2
#' @importFrom data.table ":=" setDT data.table setkey fread setnames as.data.table setcolorder melt setkeyv
#' @importFrom stats fisher.test setNames
#' @importFrom utils write.table
NULL

# Because of some issues with CRAN submission,
#(see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1") {
	utils::globalVariables(c(
	"collectionname", "collection", "filename", "size_int", "pValueLog",
	"userSet", "size", "cellType", "description", "tissue", "antibody",
	"treatment", "qValue", "oddsRatio"))
}

# This function calculates enrichment for two sets of genomic ranges intervals

# There are multiple legitimate ways to construct this kind of a test.
# It seems efficient to do what I'm doing. I take a user-set-centric approach.

# a - [support]. This could be either: the number of user regions that overlap
# at least 1 test region; or, the number of test regions that overlap
# at least 1 user region. These can be different if, for example, there are
# two test regions that overlap one user region; or vice versa. I take the
# number in the user set.
# b - [test set overlaps universe]. this could be either the number in the
# test set that overlaps at least 1 universe region, or the # in the universe
# that overlaps at least 1 test region. Using the user-set-centric approach,
# I take the # in the universe that overlap at least 1 test region.
# c - [non-hits in user set]. For this I take the size of the
# user set - support; this is the rest of the user set that did not have any
# overlaps to the test set.
# d - [size of universe - b -c -a]

# If you took a test-set centric approach, things could be slightly different.
# The user-set approach I take relies on the assumption that the user set
# regions are actually contained in the user universe. what if the user set
# contains 3 regions where the user universe is one large region? this could
# contribute 3 hits to support, but then it should also contribute 3 hits
# to b -- this could lead to errors. so, the universe must be divided on the
# divisions of the user sets at least.

#Unit Set Enrichment
#generic term "unit" refers most likely to genes, for a gene set enrichment
# analysis, but can refer to anything.
#a database set contains:
#1. a universe: the set of all possible units tested for the category of interest.
#2. the unit set: the set of all units that were positive for the test.
#3. name: an identifier that describes what this set was tested for.

# If the universe is not provided, we could consider the set of al
# ensembl genes or something like that (for gene-based)...


#a query set contains the same thing:
# 1. a query unit set: the set of units that have some property in common,
# which we wish to test for association to the database unit sets.
# 2. the tested units that were not positive.

# Both of these sets should first be restricted to the "universe"

######################################################################
# Prep functions
######################################################################

#' This function will take the user sets, overlap with the universe,
#' and redefine the user sets as the set of regions in the user
#' universe that overlap at least one region in user sets. this makes
#' for a more appropriate statistical enrichment comparison, as the user
#' sets are actually exactly the same regions found in the universe
#' otherwise, you can get some weird artifacts from the many-to-many
#' relationship between user set regions and universe regions.
#'
#' @param userSets		Regions of interest
#' @param userUniverse	Regions tested for inclusion in userSets
#' @param cores	Number of processors
#'
#' @export
#' @return	userSets redefined in terms of userUniverse
#' @example
#' R/examples/example.R
redefineUserSets = function(userSets, userUniverse, cores=1) {
	setLapplyAlias(cores)
	if(!isDisjoint(userUniverse)) {
		message("Your universe is not disjoint; try reduce() or disjoin().")
	}
	userSets =	lapplyAlias(userSets, function(x) {
		fo = findOverlaps(x, userUniverse)
		x = userUniverse[unique(subjectHits(fo))]
	})
	return(userSets)
}

#' Check universe appropriateness
#'
#' Checks to see if the universe is appropriate for the userSets
#' Anything in the userSets should be present in the universe.
#' In addition, 2 different regions in the userSets should not
#' overlap the same region in the universe
#'
#' @param userSets		Regions of interest
#' @param userUniverse	Regions tested for inclusion in userSets
#' @param cores	Number of processors
#' @param fast	Skip the (slow) test for many-to-many relationships
#'
#' @export
#' @return No return value.
#' @examples
#' data("sample_input", package="LOLA") # load userSet
#' data("sample_universe", package="LOLA") # load userUniverse
#' checkUniverseAppropriateness(userSets, userUniverse)
checkUniverseAppropriateness = function(userSets, userUniverse, cores=1, fast = FALSE) {
	message("Confirming universe appropriateness")
	userSets = listToGRangesList(userSets)
	setLapplyAlias(cores)
	userSetsLength = unlist(lapplyAlias(as.list(userSets), length))
	userSetsOlUserUniverseSum = countOverlaps(userSets, userUniverse)
	userSetsPercentInUniverseSum = userSetsOlUserUniverseSum/ userSetsLength

	if (!fast) {
		message("Checking for many-to-many relationships between sets and universe...")
		userSetsOlUserUniverseAny = countOverlapsAny(userSets, userUniverse)
		userSetsPercentInUniverseAny = userSetsOlUserUniverseAny/ userSetsLength
		message("any:", signif(userSetsPercentInUniverseAny, 6), "\n")
	} else {
		userSetsPercentInUniverseAny = 1 #skip the any test
	}

	message("sum:", signif(userSetsPercentInUniverseSum, 6), "\n")

	if (any(userSetsPercentInUniverseSum < 1)) {
		message(signif(userSetsPercentInUniverseSum, 6), "\n")
		warning(cleanws("Your user sets contain ranges that are not in your universe.
		You need to
		expand your universe. OR: your universe contains overlapping regions.
		You should reduce it. OR: your universe contains regions that overlap
		multiple regions in your user sets, You should disjoin your universe."))
#		I can check with isDisjoint)
	} else { message("PASSED") }

	if (any(userSetsPercentInUniverseAny > 1)) {
		message(signif(userSetsPercentInUniverseAny, 6), "\n")
		warning(cleanws("Your user sets contain multiple regions mapping to individual
		regions in the universe. Try redefineUserSets()"))
	} else { message("PASSED") }
}


#' If you want to test for differential enrichment within
#' your usersets, you can restrict the universe to only
#' regions that are covered in at least one of your sets.
#' This function helps you build just such a restricted
#' universe
#'
#' @param userSets The userSets you will pass to the enrichment calculation.
#' @return A restricted universe
#' @export
#' @examples

#' data("sample_input", package="LOLA") # load userSets
#' restrictedUniverse = buildRestrictedUniverse(userSets)
buildRestrictedUniverse = function(userSets) {
	disjoin(unlist(userSets))
}
