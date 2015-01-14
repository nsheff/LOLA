######################################################################
# LOADING - Functions for loading enrichment databases
######################################################################

#' Loads All Enrichment Databases.
#'
#' Helper loader functions to just load up all the data, if you want
#' to do a comprehensive analysis.
#'
#' @examples
#' loadAllEnrichmentDatabases();
#' @export
loadAllEnrichmentDatabases = function() {
	loadLocationEnrichmentDatabases();
	loadCategoryEnrichmentDatabases();
}

#' Loads all location databases. Just a helper function that calls the others.
#' use loadLocationEnrichmentMm10() for mouse.
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
	bockAnnotationMm10 <<- readRegionAnnotation(shareDir=getOption("SHARE.DATA.DIR"), bedDir="regionDB/bock_regions_mm10/");
	simpleCache("bockGRLmm10", "bockGRLmm10 = readRegionDb(bockAnnotationMm10);", cacheDir=getOption("SHARE.RCACHE.DIR"), loadEnvir=globalenv());
	bockAnnotationMm10 <<- appendAnnotations(bockAnnotationMm10, bockGRLmm10, "mm10")
}

