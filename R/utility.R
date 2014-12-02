######################################################################
# UTILITY - Small helper functions used in the package.
######################################################################

#need to get:
#replaceFileExtension();

#' setSharedDataDir
#' Sets global variable specifying the default data directory.
#' 
#' @export
setSharedDataDir = function(sharedDataDir) {
	options(SHARE.DATA.DIR=sharedDataDir); 
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

#' Just a reverser. Reverses the order of arguments and passes them untouched to countOverlapsAny -- so you can use it with lapply.
#' @export
countOverlapsAnyRev = function(subj, quer) {
	countOverlapsAny(quer, subj);
}
