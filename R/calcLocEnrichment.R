# A fork of enrichmentLocationCalc() that can use generics.
# regionDB now is a smarter object, with regions and annotation in one!
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
#' @param cores	Number of processors	
#' @param regionDB	Region DB to check for overlap, from loadRegionDB()
#' @param redefineUserSets	run redefineUserSets() on your userSets?
#'
#' @return Data.table with enrichment results
#' @export
#' @example 
#' R/examples/example.R
calcLocEnrichment = function(userSets, userUniverse, regionDB, cores=1, redefineUserSets=FALSE) {
	# Silence R CMD check Notes:
	support=d=b=userSet=pValueLog=rnkSup=rnkPV=rnkLO=NULL
	logOdds=maxRnk=meanRnk=dbSet=description=NULL
	annotationDT = regionDB$regionAnno
	testSetsGRL = regionDB$regionGRL
	annotationDT[, dbSet := 1:nrow(annotationDT)]
	setkey(annotationDT, dbSet)
	### Data sanity checks ###
	#Confirm we received GRangesList objects, convert from list if possible.
	userSets = listToGRangesList(userSets);
	testSetsGRL = listToGRangesList(testSetsGRL);
	setLapplyAlias(cores);
	
	if (any(is.null(names(testSetsGRL)))) {
		names(testSetsGRL) = 1:length(testSetsGRL);
	}

	if (redefineUserSets) { #redefine user sets in terms of universe?
		userSets =	redefineUserSets(userSets, userUniverse, cores=cores);
		userSets = listToGRangesList(userSets);
	}
	userSetsLength = unlist(lapplyAlias(as.list(userSets), length));
	
	if (! any( isDisjoint( userSets) ) ) {
		message("You have non-disjoint userSets.");
	}

	### Construct significance tests ###
	message("Calculating unit set overlaps...");
	geneSetDatabaseOverlap =lapplyAlias( as.list(userSets), countOverlapsRev, testSetsGRL);
	#geneSetDatabaseOverlap =lapplyAlias( as.list(userSets), countOverlapsAnyRev, testSetsGRL); #This is WRONG

	olmat = do.call(cbind, geneSetDatabaseOverlap); 
	#turn results into an overlap matrix. It is
	#database sets (rows) by test sets (columns), scoring the number of overlap.

	message("Calculating universe set overlaps...");
	testSetsOverlapUniverse = countOverlaps(testSetsGRL, userUniverse) #faster #returns number of items in userUniverse.
	#testSetsOverlapUniverse = countOverlapsAny(testSetsGRL, userUniverse) #returns number of items in test set
	universeLength = length(userUniverse);

	scoreTable = data.table(melt(t(olmat)))
	setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "dbSet", "support"))
	message("Calculating Fisher scores...");
	scoreTable[,c("b", "c"):=list(b=testSetsOverlapUniverse[match(dbSet, names(testSetsOverlapUniverse))]-support, c=userSetsLength-support)]
	scoreTable[,d:=universeLength-support-b-c]
	if( scoreTable[,any(b<0)] ) { #inappropriate universe.
		print(scoreTable[which(b<0),]);
		warning("Negative b entry in table. This means either: 1) Your user sets contain items outside your universe; or 2) your universe has a region that overlaps multiple user set regions, interfering with the universe set overlap calculation.");
		return(scoreTable);
		#sum(countOverlaps(testSetsGRL[[12]], userUniverse) > 0)
		#sum(countOverlaps(userUniverse, testSetsGRL[[12]]) > 0)
	}
	if( scoreTable[,any(c<0)] ) {
		warning("Negative c entry in table. Bug with userSetsLength; this should not happen.");
		return(scoreTable);
	}
	scoreTable[,c("pValueLog", "logOdds") := fisher.test(matrix(c(support,b,c,d), 2, 2), alternative='greater')[c("p.value", "estimate")], by=list(userSet,dbSet)]
	scoreTable[, pValueLog:=-log(pValueLog)]
	### Finalize and Rank results ###
	scoreTable[, rnkSup:=rank(-support, ties.method="min"), by=userSet]
	scoreTable[, rnkPV:=rank(-pValueLog, ties.method="min"), by=userSet]
	scoreTable[, rnkLO:=rank(-logOdds, ties.method="min"), by=userSet]
	scoreTable[, maxRnk:=max(c(rnkSup, rnkPV, rnkLO)), by=list(userSet,dbSet)]
	scoreTable[, meanRnk:=signif(mean(c(rnkSup, rnkPV, rnkLO)), 3), by=list(userSet,dbSet)]

	#append description column
	setkeyv(scoreTable, "dbSet")
	scoreTable = scoreTable[annotationDT]

	# limit description to 80 characters
	scoreTable[,description:=substr(description, 0, 80)]

	orderedCols = c("userSet", "dbSet", "collection", "pValueLog", "logOdds", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "cellType", "tissue", "antibody", "treatment", "dataSource", "filename", "description")
	unorderedCols = setdiff(colnames(scoreTable), orderedCols)

	setcolorder(scoreTable,  c(orderedCols, unorderedCols));
	#scoreTable[,qValue:=qvalue(pValue)$qvalue] #if you want qvalues...
	scoreTable[order(pValueLog, -meanRnk, decreasing=TRUE),]
}
