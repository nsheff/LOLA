#######################################################################
#Install and load LOLA
#######################################################################
install_github("sheffien/simpleCache") 
install_github("sheffien/LOLA") 
library(LOLA)

#######################################################################
#Populate global directory variables
#######################################################################
#Now we need to establish the connection to the shared region database.
#You must set these directory paths to match the project directory on your system! If you're on the Bock lab server, these settings will work as is.
#Nathan's Shared Repository
options(SHARE.DIR="/fhgfs/groups/lab_bock/nsheffield/share/")
source(paste0(getOption("SHARE.DIR"), "initDefault.R"))

#In the past, I would source these utilities. Now, this is deprecated as these functions are part of the LOLA package.
#utility("funcEnrichment.R")
#utility("funcGenomeLocations.R")

#This code works for two different kinds of enrichments:
#1. Locations: You give sets of genome coordinates, I look for overlaps in public datasets (right now, Encode, Cistrome, and Nathan's DNase Hypersensitivity DB).
#2. Gene Symbols: You give sets of gene symbols, I look for overlaps in public datasets (right now, mSig).

#######################################################################
#Location-style enrichments - Combined database analysis
#######################################################################
#load all databases in one go
loadLocationEnrichmentDatabases()

#your userSets is a GRangesList object OR a list object; each item in the GRangesList/list is a GRanges object with coordinates of interest (could be promoters or enhancers, for example)
#userUniverse is a GRanges object; the set of all regions tested for inclusion in your userSets.


userSets =		 	#GRangesList object
userUniverse = 		#GRanges object

#example data here:
userSets = import(paste0(getOption("SHARE.DIR"), "data/atac_example.bed"))
userUniverse =userSets


#To do a comprehensive location enrichment (using all 3 databases)
#First, you can check to see if your universe is appropriate:
checkUniverseAppropriateness(userSets, userUniverse);
locResults = locationEnrichment(userSets, userUniverse);

#View results in R:
locResults[order(pValue),][1:30,]

#Output your results:
writeCombinedEnrichment(locResults, outFolder= "locationResults", includeSplits=TRUE);




#New, Generic collection system:

regionDB = loadRegionDB(dbLocation= "~/fhgfs/share/regionDB/hg19", limit=10)
regionDB










#######################################################################
#Location-style enrichments - Individual database analysis
#######################################################################
#To do location enrichment analysis on each database individually:

#cistrome
cistromeAnnotation <<- readCistromeAnnotation(shareDir=getOption("SHARE.DATA.DIR"));
	cistromeGRL <<- simpleCache("cistromeGRL", "cistromeGRL = readCistrome(cistromeAnnotation);", cacheDir=getOption("SHARE.RCACHE.DIR"));
	cistromeAnnotation <<- appendAnnotations(cistromeAnnotation, cistromeGRL, "hg19"); #OPTIONAL
#cistromeGRL = readCistrome(cistromeAnnotation); #non-cached version
cistromeResults = enrichmentLocationCalc(userSets, userUniverse, cistromeAnnotation, cistromeGRL, dbTitle="CISTROME", mc.cores=6)

#View results sorted by p-value:
cistromeResults[order(pValue),]

#Write results to file:
writeDataTableSplitByColumn(cistromeResults, splitFactor="userSet", filePrepend="cistrome")


#encode
encodeTFBSannotation = readEncodeTFBSannotation("encodeTFBS/");
encodeGRL = simpleCache("encodeGRL", "encodeGRL = readEncodeTFBS(encodeTFBSannotation);", cacheDir=getOption("SHARE.RCACHE.DIR")); #cached version
encodeResults = enrichmentLocationCalc(userSets, userUniverse, encodeTFBSannotation, encodeGRL, dbTitle="ENCODE")

#View results sorted by p-value:
encodeResults[order(pValue),]

#Write results to file:
writeDataTableSplitByColumn(encodeResults, splitFactor="userSet", filePrepend="encode")

#dnase hypersensitivity
dhsAnnotation = readDhsAnnotation();
#dhsGRL = readDhs();#non-cached version

dhsAnnotation <<- readDhsAnnotation(shareDir=getOption("SHARE.DATA.DIR"));
dhsGRL <<- simpleCache("dhsGRL", "dhsGRL = readDhs()", cacheDir=getOption("SHARE.RCACHE.DIR"));
dhsAnnotation <<- appendAnnotations(dhsAnnotation, dhsGRL, "hg19") #optional

dhsResults = enrichmentLocationCalc(userSets, userUniverse, dhsAnnotation, dhsGRL, dbTitle="DHS", mc.cores=6)	

dhsResults[userSet==4,][order(pValue),]

writeDataTableSplitByColumn(dhsResults, splitFactor="userSet", filePrepend="dhs")

#######################################################################
#Category-style enrichments
#######################################################################
#Now your userSets and userUniverse should be list objects sets of geneSymbols.

userSets =
userUniverse = 

#mSig
loadCategoryEnrichmentDatabases()
mSigResults = enrichmentCategoryCalc(userSets, userUniverse, mSigAnnotation, mSig, dbTitle="mSig")

mSigResults[userSet==4,][order(pValue),]

writeDataTableSplitByColumn(mSigResults, splitFactor="userSet", filePrepend="mSig")

writeCombinedEnrichment()

