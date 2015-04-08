
######################################################################
# MSigDB Enrichment Functions
######################################################################

#'
#' @export
readMSigDB = function(shareDir=getOption("SHARE.DATA.DIR")) {
	#Slurp by line:
	mSig_in = readLines(paste0(shareDir, "MSigDB/msigdb.v4.0.symbols.gmt"));
	#Split on tab:
	mSig = sapply(mSig_in, stringr::str_split, "\t")
	#First entry in each line is a gene set identifier.
	names(mSig) = 1:length(mSig);
	geneSetNames = sapply(mSig, "[[", 1)
	#urls = lapply(msigdb, "[[", 2) #URLS are in spot 2, if you want 'em.

	dbLength = sapply(mSig, length);
	mSigAnnotation = data.table(clusterID=1:length(mSig), description=geneSetNames, size=dbLength)
	setkey(mSigAnnotation, "clusterID")
	#Wipe out the first two from each line
	#one is an identifier; the other is URL
	mSig = lapply(mSig, function(x) { x = x[-(1:2)]; })
	list(mSig=mSig, mSigAnnotation=mSigAnnotation)
}

#' @export
namesToNumbers = function(namedDB) {
	dbLength = sapply(namedDB, length);
	geneSetNames = names(namedDB);
	names(namedDB) = 1:length(geneSetNames);
	anno = data.table(clusterID=1:length(mSigDB), description=geneSetNames, size=dbLength)
	anno
}



#This is for Category (or annotation) enrichment calculations; for example,
#GO terms, or MSigDB
#' @export
enrichmentCategoryCalc = function(userSets, userUniverse, annotationDT, testSets, dbTitle="MSIGDB", cores=4) {
	setLapplyAlias(cores);

	if(! "list" %in% class(userSets)) {
		stop("[", dbTitle, "] userSets must be a list object.");
	}

	message("[", dbTitle, "] Calculating unit set overlaps...");
	geneSetDatabaseOverlap = lapply(userSets, function(y) { 
		sapply(testSets, function(x) { length(intersect(x, y)) } )
		} )

	olmat = do.call(cbind, geneSetDatabaseOverlap); 
	#turn results into an overlap matrix. It is
	#database sets (rows) by test sets (columns), scoring the number of overlap.

	#To make a significance test for each of these comparisons, we need to consider the total number in each dimension...

	message("[", dbTitle, "] Calculating universe set overlaps...");
	userSetsLength = sapply(userSets, length);
	testSetsOverlapUniverse = sapply(testSets, function(x) { length(intersect(x, userUniverse)) } )

	universeLength = length(userUniverse);
	scoreTable = data.table(melt(t(olmat)))

	setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "dbSet", "support")) #Sometimes this changes to X1 X2, why?
	if ( scoreTable[,class(userSet)] == "factor") {
		message("User sets were converted to factors.");
		scoreTable[,userSet:=as.character(userSet)]
	}
	message("[", dbTitle, "] Calculating Fisher scores...");
	scoreTable[,c("b", "c"):=list(b=testSetsOverlapUniverse[match(dbSet, names(testSetsOverlapUniverse))]-support, c=userSetsLength[userSet]-support)]
	scoreTable[,d:=universeLength-support-b-c]
	if( scoreTable[,any(b<0)] ) { #inappropriate universe.
		print(scoreTable[which(b<0),]);
		warning("[", dbTitle, "] Negative b entry in table. This means either: 1) Your user sets contain items outside your universe; or 2) your universe has a region that overlaps multiple user set regions, interfering with the universe set overlap calculation.");
		return(scoreTable);
		#sum(countOverlaps(testSetsGRL[[12]], userUniverse) > 0)
		#sum(countOverlaps(userUniverse, testSetsGRL[[12]]) > 0)
	}
	if( scoreTable[,any(c<0)] ) {
		stop("[", dbTitle, "] Negative c entry in table. Bug with userSetsLength; this should not happen.");
		return(scoreTable);
	}
	scoreTable[,c("pValueLog", "logOdds") := fisher.test(matrix(c(support,b,c,d), 2, 2), alternative='greater')[c("p.value", "estimate")], by=list(userSet,dbSet)]
	scoreTable[, pValueLog:=-log(pValueLog)]


	#Rank results
	scoreTable[, rnkSup:=rank(-support, ties.method="min"), by=userSet]
	scoreTable[, rnkPV:=rank(-pValueLog, ties.method="min"), by=userSet]
	scoreTable[, rnkLO:=rank(-logOdds, ties.method="min"), by=userSet]
	scoreTable[, maxRnk:=max(c(rnkSup, rnkPV, rnkLO)), by=list(userSet,dbSet)]
	scoreTable[, meanRnk:=signif(mean(c(rnkSup, rnkPV, rnkLO)), 4), by=list(userSet,dbSet)]

	#append description column
	setkeyv(scoreTable, "dbSet")
	scoreTable[annotationDT, description:=description]
	scoreTable[,db:=dbTitle]
	setcolorder(scoreTable, c("userSet", "dbSet", "description", "db", "pValueLog", "logOdds", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d"));
	scoreTable[order(pValueLog, decreasing=TRUE),]

#	#append description column
#	setkeyv(scoreTable, dbTitle)
#	scoreTable[annotationDT[, list(description=paste0(c(dbTitle, cell, treatment, antibody), collapse=" ")),by=key(annotationDT)], description:=description]
#	scoreTable
}


