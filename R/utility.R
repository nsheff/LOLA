################################################################################
# UTILITY - Small helper functions used in the package.
################################################################################
#' setSharedDataDir
#' Sets global variable specifying the default data directory.
#' 
#' @export
setSharedDataDir = function(sharedDataDir) {
	options(SHARE.DATA.DIR=sharedDataDir); 
}

countOverlapsAny = function(subj, quer, cores=1) {
	setLapplyAlias(cores)
	l = unlist(lapplyAlias(subj, function(x) { sum(overlapsAny(x, quer)) } ))
	return(l);
}
################################################################################
# INTERNAL FUNCTIONS - not exported
################################################################################

#' This will change the string in filename to have a new extension
#' @param filename	string to convert
#' @param extension	new extension
#" @returns	filename with original extension deleted, replaced by provided extension
replaceFileExtension = function(filename, extension) {
	sub("\\..*$", enforceEdgeCharacter(extension, prependChar="."), paste0(filename, "."))
}



#' Just a reverser. Reverses the order of arguments and passes them untouched to countOverlapsAny -- so you can use it with lapply.
countOverlapsAnyRev = function(subj, quer) {
	countOverlapsAny(quer, subj);
}


#' converts a list of GRanges into a GRangesList; strips all metadata.
#' @param lst	a list of GRanges objects
#' @return	a GRangesList object
listToGRangesList = function(lst) {
	if(! "GRangesList" %in% class(lst)) {
		if ("list" %in% class(lst)) {
			#strip elementMetadata
			lst = lapply(lst, function(x) { values(x) <- NULL; x; } )
			lst = GRangesList(lst);
		} else {
			warning("in listToGRangesList (funcGenomeLocations), input list must be a list object. I've taken the liberty of converting yours to a list; I hope this is OK.");
			lst = GRangesList(list(lst));
		}
	}
	return(lst);
}
