# A fork of enrichmentLocationCalc() that can use generics.
# regionDB now is a smarter object, with regions and annotation in one!
#' @export
enrichmentLocationCalcGen = function(userSets, userUniverse, regionDB, dbTitle="encode", cores=1, redefineUserSets=FALSE) {
	annotationDT = regionDB$regionAnno
	testSetsGRL = regionDB$regionGRL
	annotationDT[, dbSet := 1:nrow(annotationDT)]
	setkey(	annotationDT, dbSet)
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
	message("[", dbTitle, "] Calculating unit set overlaps...");
	geneSetDatabaseOverlap =lapplyAlias( as.list(userSets), countOverlapsRev, testSetsGRL);
	#geneSetDatabaseOverlap =lapplyAlias( as.list(userSets), countOverlapsAnyRev, testSetsGRL); #This is WRONG

	olmat = do.call(cbind, geneSetDatabaseOverlap); 
	#turn results into an overlap matrix. It is
	#database sets (rows) by test sets (columns), scoring the number of overlap.

	message("[", dbTitle, "] Calculating universe set overlaps...");
	testSetsOverlapUniverse = countOverlaps(testSetsGRL, userUniverse) #faster #returns number of items in userUniverse.
	#testSetsOverlapUniverse = countOverlapsAny(testSetsGRL, userUniverse) #returns number of items in test set
	universeLength = length(userUniverse);

	scoreTable = data.table(melt(t(olmat)))
	setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "dbSet", "support"))
	message("[", dbTitle, "] Calculating Fisher scores...");
	scoreTable[,c("b", "c"):=list(b=testSetsOverlapUniverse[match(dbSet, names(testSetsOverlapUniverse))]-support, c=userSetsLength-support)]
	scoreTable[,d:=universeLength-support-b-c]
	if( scoreTable[,any(b<0)] ) { #inappropriate universe.
		print(scoreTable[which(b<0),]);
		warning("[", dbTitle, "] Negative b entry in table. This means either: 1) Your user sets contain items outside your universe; or 2) your universe has a region that overlaps multiple user set regions, interfering with the universe set overlap calculation.");
		return(scoreTable);
		#sum(countOverlaps(testSetsGRL[[12]], userUniverse) > 0)
		#sum(countOverlaps(userUniverse, testSetsGRL[[12]]) > 0)
	}
	if( scoreTable[,any(c<0)] ) {
		warning("[", dbTitle, "] Negative c entry in table. Bug with userSetsLength; this should not happen.");
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
	scoreTable

	#append description column
	setkeyv(scoreTable, "dbSet")
	scoreTable = scoreTable[annotationDT]
	#scoreTable[,db:=dbTitle]
	setcolorder(scoreTable, c("userSet", "dbSet", "description",  "pValueLog", "logOdds", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "filename", "source", "antibody", "treatment", "collection", "size"));
	#scoreTable[,qValue:=qvalue(pValue)$qvalue] #if you want qvalues...
	scoreTable[order(pValueLog),]
}
