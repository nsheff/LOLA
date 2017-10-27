######################################################################
# ENRICHMENT - Actual workhorse enrichment calculation functions
######################################################################

#' Enrichment Calculation
#'
#' Workhorse function that calculates overlaps between userSets,
#' and then uses a fisher's exact test rank them by significance
#' of the overlap.
#'
#' @param userSets		Regions of interest
#' @param userUniverse	Regions tested for inclusion in userSets
#' @param regionDB	Region DB to check for overlap, from loadRegionDB()
#' @param minOverlap (Default:1) Minimum bases required to count an overlap
#' @param cores	Number of processors
#' @param redefineUserSets	run redefineUserSets() on your userSets?
#' @param direction    Defaults to "enrichment", but may also accept
#'     "depletion", which will swap the direction of the fisher test (use
#'     'greater' or less' value passed to the 'alternative' option of
#'     fisher.test)
#' @return Data.table with enrichment results. Rows correspond to individual
#' pairwise fisher's tests comparing a single userSet with a single databaseSet.
#' The columns in this data.table are: userSet and dbSet: index into their
#' respective input region sets. pvalueLog: -log10(pvalue) from the fisher's exact
#' result; oddsRatio: result from the fisher's exact test; support: number of
#' regions in userSet overlapping databaseSet; rnkPV, rnkOR, rnkSup: rank in this
#' table of p-value, oddsRatio, and Support respectively. The --value is the
#' negative natural log of the p-value returned from a one-sided fisher's exact
#' test. maxRnk, meanRnk: max and mean of the 3 previous ranks, providing a
#' combined ranking system. b, c, d: 3 other values completing the 2x2 contingency
#' table (with support). The remaining columns describe the dbSet for the row.
#'
#' If you have the qvalue package installed from bioconductor, runLOLA will add
#' a q-value transformation to provide FDR scores automatically.
#' @export
#' @example
#' R/examples/example.R
runLOLA = function(userSets, userUniverse, regionDB, minOverlap=1, cores=1,
redefineUserSets=FALSE, direction="enrichment") {
	# Silence R CMD check Notes:
	support=d=b=userSet=pValueLog=rnkSup=rnkPV=rnkOR=NULL
	oddsRatio=maxRnk=meanRnk=dbSet=description=NULL
	annotationDT = regionDB$regionAnno
	testSetsGRL = regionDB$regionGRL

	if (direction == "depletion") {
		fisherAlternative = "less"
	} else {
		fisherAlternative = "greater"
	}

	annotationDT[, dbSet := seq_len(nrow(annotationDT))]
	setkey(annotationDT, dbSet)
	### Data sanity checks ###
	#Confirm we received GRangesList objects, convert from list if possible.
	userSets = listToGRangesList(userSets)
	testSetsGRL = listToGRangesList(testSetsGRL)
	setLapplyAlias(cores)

	if (any(is.null(names(testSetsGRL)))) {
		names(testSetsGRL) = seq_along(testSetsGRL)
	}

	if (redefineUserSets) { #redefine user sets in terms of universe?
		userSets =	redefineUserSets(userSets, userUniverse, cores=cores)
		userSets = listToGRangesList(userSets)
	}
	userSetsLength = unlist(lapplyAlias(as.list(userSets), length))

	if (! any( isDisjoint( userSets) ) ) {
		message("You have non-disjoint userSets.")
	}

	### Construct significance tests ###
	message("Calculating unit set overlaps...")


	# Returns for each userSet, a vector of length length(testSetsGRL), with total
	# number of regions in that set overlapping anything in each testSetsGRL; this
	# is then lapplied across each userSet.


	geneSetDatabaseOverlap =
		lapplyAlias( as.list(userSets), countOverlapsRev, testSetsGRL, minoverlap=minOverlap)

	# This is WRONG:
	#geneSetDatabaseOverlap =
	#lapplyAlias( as.list(userSets), countOverlapsAnyRev, testSetsGRL)

	# This will become "support" -- the number of regions in the
	# userSet (which I implicitly assume is ALSO the number of regions
	# in the universe) that overlap anything in each database set.
	# Turn results into an overlap matrix. It is
	# dbSets (rows) by userSets (columns), counting overlap.
	olmat = do.call(cbind, geneSetDatabaseOverlap)

	message("Calculating universe set overlaps...")
	# Now for each test set, how many items *in the universe* does
	# it overlap? This will go into the calculation for c

	#faster. Returns number of items in userUniverse.
	testSetsOverlapUniverse = countOverlaps(testSetsGRL, userUniverse,
		minoverlap=minOverlap)
	# Returns number of items in test set (not used:)
	#testSetsOverlapUniverse = countOverlapsAny(testSetsGRL, userUniverse)
	# Total size of the universe
	universeLength = length(userUniverse)

	# To build the fisher matrix, support is 'a'

	scoreTable = data.table(reshape2::melt(t(olmat), variable.factor=FALSE))

	setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "dbSet", "support"))

	# reshape2 has an annoying habit of converting strings into factors, which
	# is undesirable. If the userSets are named with strings, make sure they stay
	# character. Integers are already handled appropriately.

	if ("factor" %in% class(scoreTable[, userSet])) {
		scoreTable$userSet = as.character(scoreTable$userSet)
	}

	message("Calculating Fisher scores...")
	# b = the # of items *in the universe* that overlap each dbSet,
	# less the support; This is the number of items in the universe
	# that are in the dbSet ONLY (not in userSet)
	# c = the size of userSet, less the support; This is the opposite:
	# Items in the userSet ONLY (not in the dbSet)

	scoreTable[,c("b", "c"):=list(b=testSetsOverlapUniverse[match(dbSet,
	names(testSetsOverlapUniverse))]-support, c=userSetsLength-support)]

	# This is the regions in the universe, but not in dbSet nor userSet.
	scoreTable[,d:=universeLength-support-b-c]
	if( scoreTable[,any(b<0)] ) { # Inappropriate universe.
		warning(cleanws("Negative b entry in table. This means either: 1) Your user sets
		contain items outside your universe; or 2) your universe has a region that
		overlaps multiple user set regions, interfering with the universe set overlap
		calculation."))

		return(scoreTable)
	}
	if( scoreTable[,any(c<0)] ) {
		warning("Negative c entry in table. Bug with userSetsLength; this should not happen.")
		return(scoreTable)
	}


	scoreTable[,c("pValueLog", "oddsRatio") :=

	fisher.test(matrix(c(support,b,c,d), 2, 2), alternative=fisherAlternative)[c("p.value",
	"estimate")], by=list(userSet,dbSet)]

	# Include qvalue if package exists.
	if (requireNamespace("qvalue", quietly=TRUE)) {
		# Wrap in try block since this is not vital.
		# if you want qvalues...
		tryCatch( {
			scoreTable[,qValue:=qvalue::qvalue(pValueLog)$qvalue]
		}, error = function(e) { warning("Problem in FDR calculation with qvalue.") })
	} else {
		# Another possibility for the future:
		# scoreTable[,qValue:=qValues = pmin(pValues*length(pValues),1)]
	}
	scoreTable[, pValueLog:=-log10(pValueLog)]
	### Finalize and Rank results ###
	scoreTable[, rnkSup:=rank(-support, ties.method="min"), by=userSet]
	scoreTable[, rnkPV:=rank(-pValueLog, ties.method="min"), by=userSet]
	scoreTable[, rnkOR:=rank(-oddsRatio, ties.method="min"), by=userSet]
	scoreTable[, maxRnk:=max(c(rnkSup, rnkPV, rnkOR)), by=list(userSet,dbSet)]
	scoreTable[, meanRnk:=signif(mean(c(rnkSup, rnkPV, rnkOR)), 3), by=list(userSet,dbSet)]

	# Append description column
	setkeyv(scoreTable, "dbSet")
	scoreTable = scoreTable[annotationDT]

	# limit description to 80 characters
	scoreTable[,description:=substr(description, 0, 80)]

	orderedCols = c("userSet", "dbSet", "collection", "pValueLog", "oddsRatio",
"support", "rnkPV", "rnkOR", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d",
"description", "cellType", "tissue", "antibody", "treatment", "dataSource", "filename")
	unorderedCols = setdiff(colnames(scoreTable), orderedCols)

	setcolorder(scoreTable, c(orderedCols, unorderedCols))

	scoreTable[order(pValueLog, -meanRnk, decreasing=TRUE),]
}
