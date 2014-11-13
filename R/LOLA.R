#PACKAGE DOCUMENTATION
#' Provides functions for genome location overlap analysis. Further documentation forthcoming.
#'
#' LOLA!
#' 
#' @references \url{http://github.com/sheffield}
## @import simpleCache
#' @docType package
#' @name LOLA
#' @author Nathan Sheffield
NULL

#TODO:
#move all package documentation to roxygen format.


######################################################################
# Enrichment Utility Functions
######################################################################
# These are helper functions for calculating enrichment of
# Gene sets or location clusters.
# By Nathan Sheffield, CeMM, 2014
library(reshape2)
library(data.table)

# Helper loader functions to just load up all the data, if you want
# to do a comprehensive analysis.
#' @export
loadAllEnrichmentDatabases = function() {
	loadLocationEnrichmentDatabases();
	loadCategoryEnrichmentDatabases();
}
#' @export
loadLocationEnrichmentDatabases = function() {
	#encode
	encodeTFBSannotation <<- readEncodeTFBSannotationHg19(shareDir=getOption("SHARE.DATA.DIR"));
	simpleCache("encodeGRL", "encodeGRL = readEncodeTFBS(encodeTFBSannotation);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	encodeTFBSannotation = appendAnnotations(encodeTFBSannotation, encodeGRL, "hg19")
	encodeTFBSannotation <<- appendAnnotations(encodeTFBSannotation, encodeGRL, "hg19")
	#cistrome
	cistromeAnnotation <<- readCistromeAnnotation(shareDir=getOption("SHARE.DATA.DIR"));
	simpleCache("cistromeGRL", "cistromeGRL = readCistrome(cistromeAnnotation);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	cistromeAnnotation <<- appendAnnotations(cistromeAnnotation, cistromeGRL, "hg19")
	#dnase hypersensitivity
	dhsAnnotation <<- readDhsAnnotation(shareDir=getOption("SHARE.DATA.DIR"));
	simpleCache("dhsGRL", "dhsGRL = readDhs()", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	dhsAnnotation <<- appendAnnotations(dhsAnnotation, dhsGRL, "hg19")

	message("Loaded databases: cistrome, encode, dhs.");
}

#' @export
loadCategoryEnrichmentDatabases = function() {
	#msigdb
	simpleCache("mSigList", "mSigList = readMSigDB(SHARE.RDATA.DIR)", cacheDir=getOption("SHARE.RCACHE.DIR")); 
	mSig <<- mSigList$mSig
	mSigAnnotation <<- mSigList$mSigAnnotation
	message("Loaded databases: mSig.");
}

#' @export
loadLocationEnrichmentMm9 = function() {
	encodeTFBSannotationMm9 <<- readEncodeTFBSannotationMm9(shareDir=getOption("SHARE.DATA.DIR"));
	simpleCache("encodeGRLmm9", "encodeGRLmm9 = readEncodeTFBS(encodeTFBSannotationMm9);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	cistromeAnnotationMm9 <<- readCistromeAnnotation(shareDir=getOption("SHARE.DATA.DIR"), restrictToSpecies="Mouse");
	simpleCache("cistromeGRLmm9", "cistromeGRLmm9 = readCistrome(cistromeAnnotationMm9);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
}

#' @export
loadLocationEnrichmentMm10 = function() {
	encodeTFBSannotationMm10 <<- readEncodeTFBSannotationMm10(shareDir=getOption("SHARE.DATA.DIR"));
	simpleCache("encodeGRLmm10", "encodeGRLmm10 = readEncodeTFBS(encodeTFBSannotationMm10);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	cistromeAnnotationMm10 <<- readCistromeAnnotation(shareDir=getOption("SHARE.DATA.DIR"), restrictToSpecies="Mouse");
	simpleCache("cistromeGRLmm10", "cistromeGRLmm10 = readCistrome(cistromeAnnotationMm10);", cacheDir=getOption("SHARE.RCACHE.DIR"));
	bockAnnotationMm10 <<- readRegionAnnotation(shareDir=getOption("SHARE.DATA.DIR"), bedDir="regionDB/bock_regions_mm10/", loadEnvir=globalenv());
	simpleCache("bockGRLmm10", "bockGRLmm10 = readRegionDb(bockAnnotationMm10);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	bockAnnotationMm10 <<- appendAnnotations(bockAnnotationMm10, bockGRLmm10, "mm10")
}



#helper functions
#' @export
readEncodeTFBSannotationMm9 = function(shareDir=getOption("SHARE.DATA.DIR")) {
	readEncodeTFBSannotation(encodeTFBSdir = "encodeTFBSmm9/");
}

#' @export
readEncodeTFBSannotationMm10 = function(shareDir=getOption("SHARE.DATA.DIR")) {
	readEncodeTFBSannotation(encodeTFBSdir = "encodeTFBSmm10/");
}

#' @export
readEncodeTFBSannotationHg19 = function(shareDir=getOption("SHARE.DATA.DIR")) {
	readEncodeTFBSannotation(encodeTFBSdir = "encodeTFBS/");
}


######################################################################
# ENCODE Transcription Factor Binding Site Functions
######################################################################
#These functions require a connection to SHARE.DATA.DIR, a directory
#where shared information downloaded from ENCODE is stored.

#' @export
readEncodeTFBSannotation = function(encodeTFBSdir, shareDir=getOption("SHARE.DATA.DIR")) {
	library(data.table)
	library(stringr)
	library(rtracklayer)

	#Load the files.txt file from ENCODE that annotates the
	#ChIP-seq experiments
	encodeTFBSannotation = fread(paste0(shareDir, encodeTFBSdir, "files.txt"), header=FALSE)

	#Parse the annotation into data.table format
	setnames(encodeTFBSannotation, "V1", "filename")
	ct = str_match(encodeTFBSannotation$V2, "cell=(.*?);")[, 2]
	tr = str_match(encodeTFBSannotation$V2, "treatment=(.*?);")[, 2]
	ab = str_match(encodeTFBSannotation$V2, "antibody=(.*?);")[, 2]
	encodeTFBSannotation[, V2:=NULL]
	encodeTFBSannotation[, cell:=ct]
	encodeTFBSannotation[, treatment:=tr]
	encodeTFBSannotation[, antibody:=ab]
	encodeTFBSannotation[, filename:=sub(".gz", "", filename)]
	encodeTFBSannotation[, filename:=paste0(encodeTFBSdir, "narrowPeak/", filename)]
	encodeTFBSannotation[, expID:=1:nrow(encodeTFBSannotation)] #set index variable
	setkey(encodeTFBSannotation, "expID");
	encodeTFBSannotation = getEncodeGroupSizes(encodeTFBSannotation, encodeTFBSdir);
	return(encodeTFBSannotation);
}

#' @export
getEncodeGroupSizes = function(encodeTFBSannotation, encodeTFBSdir, shareDir=getOption("SHARE.DATA.DIR")) {
	if (file.exists(paste0(shareDir, encodeTFBSdir, "groupSizes.txt"))) {
		groupSizes = fread(paste0(shareDir, encodeTFBSdir, "groupSizes.txt"))
		encodeTFBSannotation[groupSizes, size:=size]
	} else {
		encodeTFBSannotation[,size:=countFileLines(paste0(shareDir, filename)), by=expID]
		message("Recalculating and caching group sizes...");
		write.table(encodeTFBSannotation[,list(expID,size)], file=paste0(shareDir, encodeTFBSdir, "groupSizes.txt"), quote=FALSE, row.names=FALSE)
	}
	encodeTFBSannotation
}

#Load ENCODE TFBS ChIP-seq experiments into a GRangesList object.
#' @export
readEncodeTFBS = function(encodeTFBSannotation, shareDir=getOption("SHARE.DATA.DIR")) {
	grl = GRangesList()
	for (i in 1:nrow(encodeTFBSannotation)) {
		message(encodeTFBSannotation[i,]);
		filename = paste0(shareDir, "",encodeTFBSannotation$filename[[i]]);
		if (file.exists(filename)) {
			DT = fread(paste0(shareDir, "",encodeTFBSannotation$filename[[i]]))
			tfbsgr = dtToGr(DT, "V1", "V2", "V3", NULL, NULL);
			grl[[i]] = tfbsgr;
		} else {
			message("Skipping (file not found):", filename);
			grl[[i]] = GRanges();
		}
	}
	return(grl);
}

######################################################################
# CISTROME Transcription Factor Binding Site Functions
######################################################################

#Load up annotation table for encode TFBS experiments.
#' @export
readCistromeAnnotation = function(shareDir=getOption("SHARE.DATA.DIR"), cistromeDir="cistrome/", restrictToSpecies="Human") {
	library(data.table)
	library(stringr)
	library(rtracklayer)
	message("CISTROME: restricting dataset to ", restrictToSpecies);
	#Load the files.txt file from ENCODE that annotates the
	#ChIP-seq experiments
	dataCistrome = fread(paste0(shareDir, cistromeDir, "annotationsCistrome.txt"), header=TRUE)
	#Parse the annotation into data.table format
	setnames(dataCistrome, c("species", "cell", "tissue", "antibody", "treatment", "filename"))
	dataCistrome[, filename:=paste0(cistromeDir, "CistromeBeds/", filename)]

	dataEpigenome = fread(paste0(shareDir, cistromeDir, "annotationsEpigenome.txt"), header=TRUE)
	setnames(dataEpigenome, c("species", "cell", "tissue", "antibody", "treatment", "filename"))
	dataEpigenome[, filename:=paste0(cistromeDir, "EpigenomeBeds/", filename)]
	cistromeAnnotation = rbind(dataCistrome, dataEpigenome)
	
	cistromeAnnotation[,expID:=1:nrow(cistromeAnnotation)] #set index variable
	setkey(cistromeAnnotation, "expID");

	cistromeAnnotation = getCistromeGroupSizes(cistromeAnnotation, cistromeDir)
	cistromeAnnotation = cistromeAnnotation[species %in% restrictToSpecies,]
	return(cistromeAnnotation);
}

#' @export
getCistromeGroupSizes = function(cistromeAnnotation, cistromeDir="cistrome/", shareDir=getOption("SHARE.DATA.DIR")) {
	if (file.exists(paste0(shareDir, cistromeDir, "groupSizes.txt"))) {
		groupSizes = fread(paste0(shareDir, cistromeDir, "groupSizes.txt"))
		cistromeAnnotation[groupSizes, size:=size]
	} else {
		cistromeAnnotation[,size:=countFileLines(paste0(shareDir, filename, ".notrack")), by=expID]
		write.table(cistromeAnnotation[,list(expID,size)], file=paste0(shareDir, cistromeDir, "groupSizes.txt"), quote=FALSE, row.names=FALSE)
	}
	cistromeAnnotation
}

#' @export
readCistrome = function(cistromeAnnotation, shareDir=getOption("SHARE.DATA.DIR")) {
	grl = GRangesList()
	for (i in 1:nrow(cistromeAnnotation)) {
		message(cistromeAnnotation[i,]);
		filename = paste0(shareDir, cistromeAnnotation$filename[[i]]);
		if (file.exists(filename)) {
			DT = fread(paste0(filename, ".notrack"))
			tfbsgr = dtToGr(DT, "V1", "V2", "V3", NULL, NULL);
			grl[[i]] = tfbsgr;
		} else {
			message("Skipping (file not found):", filename);
			grl[[i]] = GRanges();
		}
	}
	return(grl);
}

######################################################################
# DHS Enrichment Functions
######################################################################

#need to add this into the annotation matrix:
#fread(paste0(getOption("SHARE.DATA.DIR"), "DNase/TableS05-overlapSummary.txt"))

#' @export
readDhsAnnotation = function(shareDir=getOption("SHARE.DATA.DIR")) {
	dhsDefault = data.table(clusterID=1:2500);
	setkey(dhsDefault, "clusterID")
	dhsAnno = fread(paste0(shareDir, "DNase/TableS04-cluster-to-openCellTypes.txt"))
	dhsAnnoManual = fread(paste0(shareDir, "DNase/clusterLabels.txt"), header=TRUE)
	dhsAnnotation = dhsAnno[,list(cell=paste0(unique(openTissue), collapse=";"), treatment="", antibody=""), by=clusterID]
	setkey(dhsAnnoManual, "clusterID")
	setkey(dhsAnnotation, "clusterID")
	dhsMerged = merge(dhsDefault, dhsAnnotation, all=TRUE)
	dhsMerged[is.na(cell), cell:="Weak"]
	dhsMerged = merge(dhsMerged, dhsAnnoManual, all=TRUE)
	dhsMerged[!is.na(label), cell:=label]
	dhsAnnotation = dhsMerged
	dhsSizes = fread(paste0(shareDir, "DNase/clusterSizes.txt"))
	setkey(dhsSizes, "refined_cluster");

#	dhsAnnotation[dhsSizes, size:=N]
	dhsSizes[dhsAnnotation,]
	dhsAnnotation = dhsAnnotation[dhsSizes,]

	setnames(dhsAnnotation, "N", "size")
	setnames(dhsAnnotation, "clusterID", "expID")
	return(dhsAnnotation);
}

#' @export
readDhs = function(shareDir=getOption("SHARE.DATA.DIR")) {
	message("Loading DNase database...");
	dhsClust = fread(paste0(shareDir, "DNase/TableS03-dhs-to-cluster.txt"))
	dhsClustgr = dtToGr(dhsClust, "chr", "start", "stop", NULL, NULL)
	#groupSizes = dhsClust[,.N, by=refined_cluster]
	#write.table(groupSizes, file="groupSizes.txt", quote=FALSE, row.names=FALSE)
	dhsClustList = split(dhsClustgr, dhsClust$refined_cluster)
	dhsClustList
}

######################################################################
# Generic Region DB Loading  Functions
######################################################################
#Used right now for BockDB (Mouse) 
#For a folder of bed files that has no annotation file;
#Example use:
# customAnno = readRegionAnnotation(bedDir="regionDB/diffMeth/")
# customGRL = readRegionDb(customAnno)

#' @export
readRegionAnnotation = function(shareDir=getOption("SHARE.DATA.DIR"), bedDir="regionDB/bock_regions_mm10/") { 
	DT = data.table(filename=list.files(paste0(shareDir,bedDir), "*.bed"))
	DT[, cell:=replaceFileExtension(filename, "")]
	DT[, filename:=paste0(bedDir, filename)]
	DT[, treatment:=""]
	DT[, antibody:=""]
	DT[,expID:=1:nrow(DT)]
	setkey(DT, "expID");
	DT = getRegionGroupSizes(DT, bedDir);
	DT
}

#' @export
getRegionGroupSizes = function(DT, bedDir, shareDir=getOption("SHARE.DATA.DIR")) {
	if (file.exists(paste0(shareDir, bedDir, "groupSizes.txt"))) {
		groupSizes = fread(paste0(shareDir, bedDir, "groupSizes.txt"))
		DT[groupSizes, size:=size]
	} else {
		DT[,size:=countFileLines(paste0(shareDir, filename)), by=expID]
		write.table(DT[,list(expID,size)], file=paste0(shareDir, bedDir, "groupSizes.txt"), quote=FALSE, row.names=FALSE)
	}
	DT
}

#' @export
readRegionDb = function(genericAnnotation, shareDir=getOption("SHARE.DATA.DIR")) {
	grl = GRangesList()
	for (i in 1:nrow(genericAnnotation)) {
		message(genericAnnotation[i,]);
		filename = paste0(shareDir, genericAnnotation$filename[[i]]);
		if (file.exists(filename)) {
			DT = fread(paste0(filename))
			tfbsgr = dtToGr(DT, colnames(DT)[1], colnames(DT)[2], colnames(DT)[3], NULL, NULL);
			grl[[i]] = tfbsgr;
		} else {
			message("Skipping (file not found):", filename);
			grl[[i]] = GRanges();
		}
	}
	return(grl);
}

######################################################################
# For any DB, calculate cpg coverage.
######################################################################
#This is a fantastic example of a function that caches a result the right way (april 29, 2014). I'm very proud.
#keep in mind: the name of grl is KEY!! it has to be constant.
#grl can be passed in either as a GRL (GRangesList) object, or as a character vector name of a GRL object; then I just convert in this function to grl/grlName depending on whatever you passed in. This is what lets me use the deparse/substitute method to get the name from the appendAnnotation function.
#this function should just be called with appendAnnotation, (it caches things). For a random GRL, use getCpGPercent directly in the dna methylation module. that's where the actual work is being done.
#' @export
calcCpGPercentForDb = function(grl, genomeBuild, cacheDir=getOption("SHARE.RCACHE.DIR")) {
	if(is.character(grl)) { 
		grlName = grl;
		grl = get(grlName);
	} else {
		grlName = deparse(substitute(grl));
	}
	var = paste0(grlName, "_", genomeBuild, "_cpg");
	tic();simpleCache(var, "getCpGPercent(grl, genomeBuild=genomeBuild)", buildEnvir=list(grl=grl, genomeBuild=genomeBuild), cacheDir=cacheDir);toc();
	return(get(var));
}
#anno=bockAnnotationMm10
#grl=bockGRLmm10
#' @export
appendAnnotations = function(anno, grl, genomeBuild) {
	var = deparse(substitute(grl))
	anno[,cpg:=calcCpGPercentForDb(var, genomeBuild)]
	return(anno);
}
#version for huge grls, splits them.
#appendAnnotations = function(anno, grl, genomeBuild) {
#	var = deparse(substitute(grl))
#	if ( sum(as.numeric(sum(width(encodeGRL))))

#	if ( anno[,sum(size)]  > 5e6 ) { #separate out huge ones.
#		nrow(anno)
#		splitSize = ceiling( nrow(anno) / (anno[,sum(size)] %/% 5e6) )
#		s = 1;
#		end = s + splitSize
#		i=0;
#		while (s <= nrow(anno)) {
#			i = i+1;
#			message("Split ", i);
#			currentSplit = seq(from = s, to = end);
#			print(currentSplit);
#			s = end+1;
#			end = min(s+splitSize, nrow(anno));
#			var_update = paste0(var, "_split", i);
#			message("Size: ", anno[currentSplit,sum(size)])
#			message(var_update);
#			assign(var_update, grl[currentSplit], env=.GlobalEnv)
#			anno[currentSplit,cpg:=calcCpGPercentForDb(var_update, genomeBuild)]
#		}
#	} else {
#		anno[,cpg:=calcCpGPercentForDb(var, genomeBuild)]
#	}
#	return(anno);
#}


######################################################################
# MSigDB Enrichment Functions
######################################################################

#' @export
readMSigDB = function(shareDir=getOption("SHARE.DATA.DIR")) {
	library(stringr)
	#Slurp by line:
	mSig_in = readLines(paste0(shareDir, "MSigDB/msigdb.v4.0.symbols.gmt"));
	#Split on tab:
	mSig = sapply(mSig_in, str_split, "\t")
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


######################################################################
# HELPER FUNCTIONS 
######################################################################

#two helper functions
#' @export
countOverlapsAny = function(subj, quer, cores=1) {
	setLapplyAlias(cores)
	l = unlist(lapplyAlias(subj, function(x) { sum(overlapsAny(x, quer)) } ))
	return(l);
}

#' @export
countOverlapsAnyRev = function(subj, quer) {
	countOverlapsAny(quer, subj);
}


######################################################################
# ENRICHMENT - Actual workhorse enrichment calculation functions
######################################################################

# convenience function for doing all the location enrichment functions
# in one fell swoop.
#' @export
locationEnrichment = function(userSets, userUniverse, checkUniverse=FALSE, cores=4, redefineUserSets=TRUE) {
	if (checkUniverse & !redefineUserSets) {
		checkUniverseAppropriateness(userSets, userUniverse, cores);
	}
	cistromeResults = enrichmentLocationCalc(userSets, userUniverse, cistromeAnnotation, cistromeGRL, dbTitle="CISTROME", cores=cores, redefineUserSets=redefineUserSets)
	encodeResults = enrichmentLocationCalc(userSets, userUniverse, encodeTFBSannotation, encodeGRL, dbTitle="ENCODE", cores=cores, redefineUserSets=redefineUserSets)
	dhsResults = enrichmentLocationCalc(userSets, userUniverse, dhsAnnotation, dhsGRL, dbTitle="DHS", cores=cores, redefineUserSets=redefineUserSets)
	combinedResults = rbind(cistromeResults, encodeResults, dhsResults)
	return(combinedResults[order(pValueLog, decreasing=TRUE),]);
}


#this function will take the user sets, overlap with the universe, and redefine the user sets as the set of regions in the user universe that overlap at least one region in user sets. this makes for a more appropriate statistical enrichment comparison, as the user sets are actually exactly the same regions found in the universe; otherwise, you can get some weird artifacts from the many-to-many relationship between user set regions and universe regions.
#' @export
redefineUserSets = function(userSets, userUniverse, cores=1) {
	setLapplyAlias(cores);
	if(!isDisjoint(userUniverse)) {
		message("Your universe is not disjoint; try reduce() or disjoin().");
	}
	userSets =	lapplyAlias(userSets, function(x) { fo = findOverlaps(x, userUniverse); x = userUniverse[unique(subjectHits(fo))]; } )
	return(userSets);
}
0

#' @export
checkUniverseAppropriateness = function(userSets, userUniverse, cores=1, fast = FALSE) {
	message("Confirming universe appropriateness");
	userSets = listToGRangesList(userSets);
	setLapplyAlias(cores)
	userSetsLength = unlist(lapplyAlias(as.list(userSets), length));
	userSetsOlUserUniverseSum = countOverlaps(userSets, userUniverse); 
	userSetsPercentInUniverseSum = userSetsOlUserUniverseSum/ userSetsLength;

	if (!fast) {
		message("Checking for many-to-many relationships between sets and universe...");
		userSetsOlUserUniverseAny = countOverlapsAny(userSets, userUniverse); 
		userSetsPercentInUniverseAny = userSetsOlUserUniverseAny/ userSetsLength;
		cat("any:", signif(userSetsPercentInUniverseAny, 6), "\n");
	} else {
		userSetsPercentInUniverseAny = 1; #skip the any test;
	}

	cat("sum:", signif(userSetsPercentInUniverseSum, 6), "\n");

	if (any(userSetsPercentInUniverseSum < 1)) {
		cat(signif(userSetsPercentInUniverseSum, 6), "\n");
		warning("Your user sets contain ranges that are not in your universe. You need to expand your universe. OR: your universe contains overlapping regions. You should reduce it.");
#		I can check with isDisjoint)
	}

	if (any(userSetsPercentInUniverseAny > 1)) {
		cat(signif(userSetsPercentInUniverseAny, 6), "\n");
		warning("Your user sets contain multiple regions mapping to individual regions in the universe. Try redefineUserSets()");
	}

}

#This function calculates enrichment for two sets of genomic ranges intervals

#There are multiple legitimate ways to construct this kind of a test. It seems most efficient to do what I'm doing. I take a user-set-centric approach.

#a - [support]. This could be either: the number of user regions that overlap at least 1 test region; or, the number of test regions that overlap at least 1 user region. These can be different if, for example, there are two test regions that overlap one user region; or vice versa. I take the # in the user set.
#b - [test set overlaps universe]. this could be either the # in the test set that overlaps at least 1 universe region, or the # in the universe that overlaps at least 1 test region. Using the user-set-centric approach, I take the # in the universe that overlap at least 1 test region.
#c - [non-hits in user set]. For this I take the size of the user set - support; this is the rest of the user set that did not have any overlaps to the test set.
#d - [size of universe - b -c -a]

#if you took a test-set centric approach, things could be slightly different. The user-set approach I take relies on the assumption that the user set regions are actually contained in the user universe. what if the user set contains 3 regions where the user universe is one large region? this could contribute 3 hits to support, but then it should also contribute 3 hits to b -- this could lead to errors. so, the universe must be divided on the divisions of the user sets at least.

#' @export
enrichmentLocationCalc = function(userSets, userUniverse, annotationDT, testSetsGRL, dbTitle="encode", cores=1, redefineUserSets=FALSE) {
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
	scoreTable[annotationDT[, list(description=paste0(c(cell, treatment, antibody), collapse=" ")),by=key(annotationDT)], description:=description]
	scoreTable[,db:=dbTitle]
	setcolorder(scoreTable, c("userSet", "dbSet", "description", "db", "pValueLog", "logOdds", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d"));
	#scoreTable[,qValue:=qvalue(pValue)$qvalue] #if you want qvalues...
	scoreTable[order(pValueLog),]
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







######################################################################
# POST-ENRICHMENT - Functions for processing enrichment results
######################################################################

#this function, given a single row from an enrichment table calculation, it will find the set of overlaps between the user set and the test set. You can then use these, for example, to get fasta sequences for those regions.
#' @export
extractEnrichmentOverlaps = function(locResult, userSets) {
	print(locResult);
	dbGRL = switch(locResult[, db], 
		"ENCODE" = encodeGRL,
		"CISTROME" = cistromeGRL,
		"DHS" = dhsGRL)
	message(length(dbGRL));
	userSet = userSets[[as.character(locResult[,userSet])]]
	dbSet = locResult[1,dbSet]
	if (any(!is.null(names(dbGRL)))) { #convert dbset named into number
		dbSet = match(dbSet, names(dbGRL)); 
	}
	userSet[queryHits(findOverlaps(userSet,dbGRL[[dbSet]]))]
}

#efficiently split a data.table
#' @export
splitDataTable = function(DT, splitFactor) {
	split(1:nrow(DT), DT[, get(splitFactor)])
}



#Given a data table and a factor variable to split on,
#efficiently divides the table and then writes the different splits
#to separate files, named with filePrepend and numbered according
#to split.
#@returns number of splits written
#' @export
writeDataTableSplitByColumn = function(DT, splitFactor, filePrepend="", orderColumn=NULL) {
	saveScipenSetting = getOption("scipen"); 
	options(scipen = 4); #use scientific notation for pvalues.
	if (is.null(orderColumn)) {
		orderColumn = colnames(DT)[1];	#default order by first col.
	}
	length(
lapply( split(1:nrow(DT), DT[, get(splitFactor)]), 
		function(x) {
			fileName = paste0(filePrepend, DT[x,get(splitFactor)][1], ".txt");
			if (file.exists(fileName)) {
				message("Overwriting ", fileName , "...");
			} else {
				message(fileName);
			}
			write.table(DT[x,][order(get(orderColumn)),], file=fileName, quote=FALSE, row.names=FALSE, sep="\t"); 
		}
	)
	);
	options(scipen = saveScipenSetting);
}


#Function for writing output all at once: combinedResults is an table 
#generated by "locationEnrichment()".# or by rbinding category and location results.
#Writes all enrichments to a single file, and also spits out the same data
#divided into groups based on userSets, and Databases, just for convenience.
#disable this with an option.
#' @export
writeCombinedEnrichment = function(combinedResults, outFolder=NULL, includeSplits=TRUE) {
	if (outFolder == "" | is.null(outFolder)) {
		outFolder = "";
	} else if (substr(outFolder, nchar(outFolder), nchar(outFolder)) != "/") {
		outFolder = paste0(outFolder, "/");
	}

	dir.create(outFolder, showWarnings=FALSE);
	if (includeSplits) {
		if (combinedResults[,length(unique(userSet))] > 1) {
			writeDataTableSplitByColumn(combinedResults[order(pValueLog,decreasing=TRUE),], splitFactor="userSet", filePrepend=paste0(outFolder, "userSet_"));
		}
		if (combinedResults[,length(unique(db))] > 1) {
			writeDataTableSplitByColumn(combinedResults[order(pValueLog,decreasing=TRUE),], splitFactor="db", filePrepend=paste0(outFolder, "db_"));
		}
	}
	if (file.exists(paste0(outFolder, "allEnrichments.txt")))
		message("Overwriting ", paste0(outFolder, "allEnrichments.txt"), "...");
	write.table(combinedResults[order(pValueLog,decreasing=TRUE),], file=paste0(outFolder, "allEnrichments.txt"), row.names=FALSE, quote=FALSE, sep="\t");
}
















# Deprecated function

getTopEnrichedHits = function(sigvals, n, annotationTable=NULL) {
	topn = order(sigvals)[1:n]
	if (is.null(annotationTable)) {
		return(	data.table(cbind(dbGeneSet=topn, pval=sigvals[topn])) )
	}
		return (cbind(category=topn, pval=sigvals[topn],annotationTable[topn]) )
}











#Unit Set Enrichment
#generic term "unit" refers most likely to genes, for a gene set enrichment analysis, but can refer to anything.
#a database set contains:
#1. a universe: the set of all possible units that were tested for the category of interest.
#2. the unit set: the set of all units that were positive for the test.
#3. name: an identifier that describes what this set was tested for.

#If the universe is not provided, we could consider the set of all ensembl genes or something like that (for gene-based)...


#a query set contains the same thing:
#1. a query unit set: the set of units that have some property in common, which we wish to test for association to the database unit sets.
#2. the tested units that were not positive.

#Both of these sets should first be restricted to the "universe"



#DEPRECATED FUNCTION, SPECIFIC TO ENCODE. 
#I RE-WROTE THIS FUNCTION TO ENABLE IT TO WORK WITH ANY DATABASE, 
#SUCH AS CISTROME.
# Using the ENCODE TFBS annotations:

#encodeTFBSannotation = readEncodeTFBSannotation(cacheDir);
#encodeGRL = simpleCache("encodeGRL", "encodeGRL = readEncodeTFBS(encodeTFBSannotation);", cacheDir=SHARE.RDATA.DIR);
#r = countOverlaps(encodeGRL, mygr)


#function that calculates encode overlap
#gRangesList should be a list of GRanges objects,
#which are, for example, the PROMOTERS of the gene sets you are interested in
#enrichmentEncodeTFBS = function(gRangesSetsList, userUniverse, mc.cores=1) {
#	encodeTFBSannotation = readEncodeTFBSannotation(cacheDir);
#	if (is.function(simpleCache)) {
#		encodeGRL = simpleCache("encodeGRL", "encodeGRL = readEncodeTFBS(encodeTFBSannotation);", cacheDir=SHARE.RDATA.DIR);
#	} else {
#		message("Loading ENCODE data...");
#		encodeGRL = readEncodeTFBS(encodeTFBSannotation);
#	}

#	if(mc.cores > 1) { #multicore for this? could take a while...
#		library(multicore)
#		lapplyAlias = mclapply;
#	} else {
#		lapplyAlias = lapply;
#	}
#	encodeGRL.length = lapplyAlias(encodeGRL, length, mc.cores=6);
#	message("Calculating unit set overlaps...");
#	geneSetDatabaseOverlap =lapplyAlias( as.list(gRangesSetsList), countOverlapsRev, encodeGRL, mc.cores=6);

#	olmat = do.call(cbind, geneSetDatabaseOverlap); 
#	#turn results into an overlap matrix. It is
#	#database sets (rows) by test sets (columns), scoring the number of overlap.

#	#To make a significance test for each of these comparisons, we need to consider the total number in each dimension...

#	message("Calculating universe set overlaps...");
#	gRangesSetsList.length = sapply(gRangesSetsList, length);
#	#univPromOverlaps =countOverlaps(userUniverse, encodeGRL);
#	encOverlaps =countOverlaps(encodeGRL,userUniverse);

#	universeLength = length(userUniverse);

#	scoreTable = data.table(melt(t(olmat)))
#	setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "expID", "support"))

#	message("Calculating Fisher scores...");
#	scoreTable[,c("b", "c", "d"):=list(b=encOverlaps[expID]-support, c=gRangesSetsList.length[userSet]-support, d=universeLength-encOverlaps[expID])]
#	scoreTable[,c("pValue", "logOdds") := fisher.test(matrix(c(support,b,c,d), 2, 2), alternative='greater')[c("p.value", "estimate")], by=list(userSet,expID)]
#	scoreTable

#	message("Ranking results...");
#	scoreTable[, rankSupport:=rank(-support, ties.method="min"), by=userSet]
#	scoreTable[, rankPValue:=rank(pValue, ties.method="min"), by=userSet]
#	scoreTable[, rankLogOdds:=rank(-logOdds, ties.method="min"), by=userSet]
#	scoreTable[, maxRank:=max(c(rankSupport, rankPValue, rankLogOdds)), by=list(userSet,expID)]
#	scoreTable[, meanRank:=signif(mean(c(rankSupport, rankPValue, rankLogOdds)), 4), by=list(userSet,expID)]

#	#append description column
#	setkey(scoreTable, "expID")
#	scoreTable[encodeTFBSannotation[, list(description=paste0(c("ENCODE", cell, treatment, antibody), collapse=" ")),by=expID], description:=description]
#	scoreTable
#}




