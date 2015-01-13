######################################################################
# Generic Region DB Loading  Functions
######################################################################
#' Helper function to annotate and load a regionDB, a folder with
#' subfolder collections of regions.
#' 
#' @param dbLocation	folder where your regionDB is stored.
#' @param filePattern	passed to list.files; you can use this to select only certain file names in your folders.
#' @export
#' @examples
#' regionDB = loadRegionDB(dbLocation= "~/fhgfs/share/regionDB/hg19")
loadRegionDB = function(dbLocation, filePattern="") {
	regionAnno = readRegionAnnotation(dbLocation, filePattern);
	regionGRL = readRegionGRL(dbLocation, regionAnno);
	return(nlist(dbLocation, regionAnno, regionGRL));
}

#' Given a folder containing region collections in subfolders, this function
#' will either read the annotation file if one exists, or create a generic
#' annotation file.
#' 
#' @export
readRegionAnnotation = function(dbLocation, filePattern = "") {
	if (is.null(shareDir)) {
		message("You must set global option SHARE.DATA.DIR with setSharedCacheDir(), or specify a shareDir parameter directly to readBeds().");
		return(NA);
	}
	#Build a data.table annotating the beds.
	#Should give collections
	collections = list.dirs(path=dbLocation, full.names=FALSE, recursive=FALSE)
	message("Found collections: ", paste(collections, collapse=", "));
	annoDT = data.table();
	annotationColNames = c("filename", "description", "source", "antibody", "treatment")

	for (collection in collections) {
		files = list.files(paste0(dbLocation,"/",collection), filePattern)
		#eliminate special annotation files
		files = files [ -grep("^0", files)]
		collectionAnnoDT = data.table(collection=collection, filename=files); #preserve new ones
		setkey(collectionAnnoDT, "filename")

		#look for index file
		indexFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "0index")
		if (file.exists(indexFile)) {
			indexDT = fread(paste0(dbLocation, collection, "0index"));
			setcolnames(indexDT, lc(colnames(indexDT)));
			indexDT = indexDT[,annotationColNames] #subset
		} else {
			indexDT = as.data.table(setNames(replicate(length(annotationColNames),character(0), simplify = F), annotationColNames));
			setkey(indexDT, "filename")
		}
		#look for size file
		sizeFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "0sizes")
		if (file.exists(sizeFile)) {
			groupSizes = fread(sizeFile)
			#setkey(groupSizes, filename);
			collectionAnnoDT[groupSizes, size:=size]
		} else {
			message("Collection: ", collection, ". Creating size file...")
			collectionAnnoDT[,size:=countFileLines(paste0(dbLocation, "/", collection, "/", filename)), by=filename]
			write.table(collectionAnnoDT[,list(filename,size)], file=sizeFile, quote=FALSE, row.names=FALSE, sep="\t")
		}
	collectionAnnoDT = indexDT[collectionAnnoDT]
	annoDT = rbind(annoDT, collectionAnnoDT);
	} #end loop through collections

	return(annoDT)
}

#' Given a region annotation object, this function either reads its sizes
#' file if one exists, or otherwise calculates sizes and creates a sizes
#' file.
#"
#' @param dbLocation	folder of regiondB
#' @param annoDT 	the result of readRegionAnnotation();
#' @export
getRegionGroupSizes = function(dbLocation, annoDT) {
	dbLocation = enforceTrailingSlash(dbLocation);
	collections = annoDT[,unique(collection)]
	for (collection in collections) {
		sizeFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "0sizes")
		if (file.exists(sizeFile)) {
			groupSizes = fread(sizeFile)
			#setkey(groupSizes, filename);
			annoDT[groupSizes, size:=size]
		} else {
			message("Collection: ", collection, ". Creating size file...")
			annoDT[,size:=countFileLines(paste0(dbLocation, "/", collection, "/", filename)), by=filename]
			write.table(annoDT[,list(filename,size)], file=sizeFile, quote=FALSE, row.names=FALSE, sep="\t")
		}
	}
	annoDT
}

#' This function takes a region annotation object and reads in the regions,
#' returning a GRangesList object of the regions.
#' 
#' @param dbLocation	folder of regiondB
#' @param annoDT	output of readRegionAnnotation().
#' @param limit	for testing purposes, you could limit the number of files read. NULL for no limit (default).
#' @export
readRegionGRL = function(dbLocation, annoDT, limit=NULL) {
	grl = GRangesList()
	dbLocation = enforceTrailingSlash(dbLocation);
	filesToRead = annoDT[,list(fullFilename=paste0(dbLocation, enforceTrailingSlash(collection), filename)), by=filename]$fullFilename

	if (is.null(limit)) {
		limit = length(filesToRead);
	}	
	for (i in 1:limit) {
		message(i, ": ", filesToRead[i]);
		filename = filesToRead[i]
		if (file.exists(filename)) {
			success = tryCatch( { 
				DT = fread(paste0(filename))
				tfbsgr = dtToGr(DT, colnames(DT)[1], colnames(DT)[2], colnames(DT)[3], NULL, NULL);
				grl[[i]] = tfbsgr;
				TRUE
			},
			error = function(e) { message(i, " ERR:", filename); return(FALSE); } )
			if (!success) { grl[[i]] = GRanges(); }
		} else {
			message("Skipping (file not found):", filename);
			grl[[i]] = GRanges();
		}
	}
	return(grl);
}


