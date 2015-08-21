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
#' @param sharedDataDir	directory where the shared data is stored.
#' @return No return value.
#' @examples
#' setSharedDataDir("project/data")
setSharedDataDir = function(sharedDataDir) {
	options(SHARE.DATA.DIR=sharedDataDir)
	return()
}


countOverlapsAny = function(subj, quer, cores=1) {
	setLapplyAlias(cores)
	l = unlist(lapplyAlias(subj, function(x) { sum(overlapsAny(x, quer)) } ))
	return(l)
}
################################################################################
# INTERNAL FUNCTIONS - not exported
################################################################################

#' This will change the string in filename to have a new extension
#' @param filename	string to convert
#' @param extension	new extension
#' @return	Filename with original extension deleted, replaced by provided
#' extension
replaceFileExtension = function(filename, extension) {
	sub("\\..*$", enforceEdgeCharacter(extension, prependChar="."),
	paste0(filename, "."))
}



#' Just a reverser. Reverses the order of arguments and passes them
#' untouched to countOverlapsAny -- so you can use it with lapply.
#'
#' @param subj Subject
#' @param quer Query
#' @return Results from countOverlaps
countOverlapsAnyRev = function(subj, quer) {
	countOverlapsAny(quer, subj)
}


#' converts a list of GRanges into a GRangesList; strips all metadata.
#' @param lst	a list of GRanges objects
#' @return	a GRangesList object
listToGRangesList = function(lst) {
	if(! "GRangesList" %in% class(lst)) {
		if ("list" %in% class(lst)) {
			#strip elementMetadata
			lst = lapply(lst, function(x) { values(x) <- NULL; x; } )
			lst = GRangesList(lst)
		} else {
			warning("in listToGRangesList (funcGenomeLocations), input list must be
			a list object. I've taken the liberty of converting yours to a list; I
			hope this is OK.")
			lst = GRangesList(list(lst))
		}
	}
	return(lst)
}


#' Wrapper of write.table that provides defaults to write a
#' simple .tsv file. Passes additional arguments to write.table
#'
#' @param ... Additional arguments passed to write.table
#' @return No return value
write.tsv = function(...) {
	write.table(..., sep="\t", row.names=FALSE, quote=FALSE)
}

######################################################################

#' Imports bed files and creates GRanges objects, using the fread()
#' function from data.table.
#'
#' @param file File name of bed file.
#' @return GRanges Object
#' @export
#' @examples
#' a = readBed(system.file("extdata", "examples/combined_regions.bed",
#' package="LOLA"))
readBed = function(file) {
	DT = fread(file)
	# bed specification says:
	# 1=chr, 2=start, 3=end, 4=name, 5=score (discarded), 6=strand.
	cn = rep(NA, 6)
	readCols = colnames(DT)
	cn[seq_len(readCols)] = readCols
	tfbsgr = dtToGr(DT, chr=cn[1], start=cn[2], end=cn[3],
		name=cn[4], strand=cn[6])
	return(tfbsgr)
}


#Two utility functions for converting data.tables into GRanges objects
#genes = dtToGR(gModels, "chr", "txStart", "txEnd", "strand", "geneId")
dtToGrInternal = function(DT, chr, start, end=NA,
	strand=NA, name=NA, metaCols=NA) {
	if (is.na(end)) {
		end = start
	}
	if (is.na(strand)) {
			gr=GRanges(
				seqnames=DT[[`chr`]],
				ranges=IRanges(start=DT[[`start`]],
				end=DT[[`end`]]),
				strand="*"
			)
	} else {
		gr=GRanges(
			seqnames=DT[[`chr`]],
			ranges=IRanges(start=DT[[`start`]],
			end=DT[[`end`]]),
			strand=DT[[`strand`]]
		)
	}
	if (! is.na(name) ) {
		names(gr) = DT[[`name`]]
	} else {
		names(gr) = seq_along(gr)
	}
	if(! is.na(metaCols)) {
		for(x in metaCols) {
			elementMetadata(gr)[[`x`]]=DT[[`x`]]
		}
	}
	gr
}

dtToGr = function(DT, chr="chr", start="start", end=NA, strand=NA, name=NA,
	splitFactor=NA, metaCols=NA) {
	if(is.na(splitFactor)) {
		return(dtToGrInternal(DT, chr, start, end, strand, name,metaCols))
	}
	if ( length(splitFactor) == 1 ) {
		if( splitFactor %in% colnames(DT) ) {
			splitFactor = DT[, get(splitFactor)]
		}
	}
	lapply(split(seq_len(nrow(DT)), splitFactor),
			function(x) {
				dtToGrInternal(DT[x,], chr, start, end, strand, name,metaCols)
			}
		)


}

# If you want to use the GenomicRanges countOverlaps function, but you want to
# do it in an lapply, that will work... but you can only do it in one direction.
# If you want to lapply on the opposite argument, you can't do it (because
# countOverlaps is not symmetric: it depends on which argument comes first).
# If you want to do an lapply, but countOverlaps with the query as the second
# argument instead of the first, you can use this function to simply reverse
# the order of the arguments.
# This is used in the enrichment calculations (originally from the EWS
# project; 2014, CeMM).
countOverlapsRev = function(query, subject) {
	return(countOverlaps(subject, query))
}



# Parses result of system "wc" wordcount to return the number of lines
# in a file into R.
countFileLines = function(filename) {
	if (!file.exists(filename)) {
		warning("File does not exist:", filename); return(0);
	}
	as.numeric(
		strsplit(sub("^\\s+", "",
		system(paste("wc -l ", filename),
		intern=TRUE)), " ")[[1]][1]
	)
}


#' Function to sample regions from a GRangesList object, in specified proportion
#'
#' @param GRL	GRangesList from which to sample
#' @param prop	vector with same length as GRL, of values between 0-1,
#' proportion of the list to select
#'
#' @return A sampled subset of original GRangesList object.
sampleGRL = function(GRL, prop) {
	sampleGRanges = function(GR, prop) {
		GR[sample(length(GR), floor(length(GR) * prop))]
	}
	mapply(sampleGRanges, GRL, prop)
}


#' To make parallel processing a possibility but not required,
#' I use an lapply alias which can point at either the base lapply
#' (for no multicore), or it can point to mclapply,
#' and set the options for the number of cores (what mclapply uses).
#'
#' @param cores	Number of cpus
#' @return No return value
#' @export
#' @examples
#' setLapplyAlias(4)
#' setLapplyAlias(1)
setLapplyAlias = function(cores) {
	if(cores > 1) { #use multicore?
	if (requireNamespace("parallel", quietly = TRUE)) {
		options(mc.cores=cores)
	} else {
		warning("You don't have package parallel installed. Setting cores to 1.")
		options(mc.cores=1); #reset cores option.
		}
	} else {
		options(mc.cores=1); #reset cores option.
	}
	return()
}

#' Function to run lapply or mclapply, depending on the option set in
#' getOption("mc.cores"), which can be set with setLapplyAlias().
#'
#' @param ... Arguments passed lapply() or mclapply()
#' @return Result from lapply 0r parallel::mclapply
#' @export
#' @examples
#' lapplyAlias(letters, paste0, ".")
lapplyAlias = function(...) {
	if(getOption("mc.cores") > 1) {
		return(parallel::mclapply(...))
	} else {
		return(lapply(...))
	}
}





#check for, and fix, trailing slash. if necessary
enforceTrailingSlash = function(folder) {
	enforceEdgeCharacter(folder, appendChar="/")
}
enforceEdgeCharacter = function(string, prependChar="", appendChar="") {
	if (string=="" | is.null(string)) {
		return(string)
	}
	if(!is.null(appendChar)) {
		if (substr(string,nchar(string), nchar(string)) != appendChar) { # +1 ?
			string = paste0(string, appendChar)
			}
	}
	if (!is.null(prependChar)) {
		if (substr(string,1,1) != prependChar) { # +1 ?
			string = paste0(prependChar, string)
		}
	}
	return(string)
}


#' Nathan's magical named list function.
#' This function is a drop-in replacement for the base list() function,
#' which automatically names your list according to the names of the
#' variables used to construct it.
#' It seemlessly handles lists with some names and others absent,
#' not overwriting specified names while naming any unnamed parameters.
#' Took me awhile to figure this out.
#'
#' @param ...	arguments passed to list()
#' @return A named list object.
nlist = function(...) {
	fcall = match.call(expand.dots=FALSE)
	l = list(...)
	if(!is.null(names(list(...)))) {
		names(l)[names(l) == ""] = fcall[[2]][names(l) == ""]
	} else {
		names(l) = fcall[[2]]
	}
	return(l)
}
