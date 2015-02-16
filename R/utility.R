################################################################################
# UTILITY FUNCTIONS
################################################################################
# These are functions copied over from my repository of utilities used
# by this package. They are repeated here simply for portability, so this
# package can be deployed on systems without access to my utilities. 
# Any changes should probably be backported to the primary functions rather 
# than in these convenience duplications.
#
# These functions should probably remain interior to the package (not exported)
#
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

######################################################################
#Two utility functions for converting data.tables into GRanges objects
#genes = dtToGR(gModels, "chr", "txStart", "txEnd", "strand", "geneId");
dtToGrInternal = function(DT, chr, start, end=NULL, str=NULL, name=NULL,metaCols=NULL) {
	if (is.null(end)) {
		end = start;
	}
	if (is.null(str)) {
			gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), str="*")
	} else {
	gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), str=DT[[`str`]])
	}
	if (! is.null(name)) {
		names(gr) = DT[[`name`]];
	} else {
		names(gr) = 1:length(gr);
	}
	if(! is.null(metaCols)) {
		for(x in metaCols) {
			elementMetadata(gr)[[`x`]]=DT[[`x`]]
		}
	}
	gr;
}

dtToGr = function(DT, chr="chr", start="start", end=NULL, strand=NULL, name=NULL, splitFactor=NULL, metaCols=NULL) {
	if(is.null(splitFactor)) {
		return(dtToGrInternal(DT, chr, start, end, strand, name,metaCols));
	}
	if ( length(splitFactor) == 1 ) { 
		if( splitFactor %in% colnames(DT) ) {
			splitFactor = DT[, get(splitFactor)];
		}
	}
	lapply(split(1:nrow(DT), splitFactor), 
			function(x) { 
				dtToGrInternal(DT[x,], chr, start, end, strand, name,metaCols)
			}
		)


}

#If you want to use the GenomicRanges countOverlaps function, but you want to do it in an lapply, that will work... but you can only do it in one direction. If you want to lapply on the opposite argument, you can't do it (because countOverlaps is not symmetric: it depends on which argument comes first). If you want to do an lapply, but countOverlaps with the query as the second argument instead of the first, you can use this function to simply reverse the order of the arguments.
#This is used in the enrichment calculations (originally from the EWS project; 2014, CeMM).
countOverlapsRev = function(query, subject) {
	return(countOverlaps(subject, query));
}



# Parses result of system "wc" wordcount to return the number of lines in a file into R.
countFileLines = function(filename) {
	if (!file.exists(filename)) { warning("File does not exist:", filename); return(0); }
	as.numeric(strsplit(system(paste("wc -l ", filename), intern=TRUE), " ")[[1]][1])
}



#To make multicore a possibility but not required, I use an lapply alias which can point at either the base lapply (for no multicore), or it can load library(multicore) and then point to mclapply, and set the options for the number of cores (which is what mclapply uses).
setLapplyAlias = function(cores) {
	if(cores > 1) { #use multicore?
		library(parallel)
		options(mc.cores=cores);
		lapplyAlias <<- mclapply;
	} else {
		lapplyAlias <<- lapply;
		options(mc.cores=1); #reset cores option.
	}
}





