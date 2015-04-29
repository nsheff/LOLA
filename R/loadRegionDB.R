######################################################################
# Generic Region DB Loading  Functions
######################################################################


#' Helper function to annotate and load a regionDB, a folder with
#' subfolder collections of regions.
#' 
#' @param dbLocation	folder where your regionDB is stored.
#' @param filePattern	passed to list.files; you can use this to select 
#'	only certain file names in your folders.
#' @param useCache	uses simpleCache to cache and load the results
#' @param limit 	You can limit the number of regions for testing.
#'	Default: NULL (no limit)
#'
#' @return regionDB list containing database location, region and 
#' collection annotations, and regions GRangesList
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionDB = loadRegionDB(dbLocation=dbPath)
loadRegionDB = function(dbLocation, filePattern="", useCache=TRUE, limit=NULL) {
	regionAnno = readRegionSetAnnotation(dbLocation, filePattern);
	collectionAnno = readCollectionAnnotation(dbLocation);
	regionGRL = readRegionGRL(dbLocation, regionAnno, useCache, limit=limit);
	return(nlist(dbLocation, regionAnno, collectionAnno, regionGRL));
}

#' Read collection annotation
#'
#' @param dbLocation	Location of the database
#'
#' @return Collection annotation data.table
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' collectionAnno = readCollectionAnnotation(dbLocation=dbPath)
readCollectionAnnotation = function(dbLocation) {
	annoDT = data.table();
	collections = list.dirs(path=dbLocation, full.names=FALSE, recursive=FALSE)
	message("Reading collection annotations: ", paste(collections, collapse=", "));
	collectionColNames = c("collector", "date", "source", "description")
	collectionsDT = data.table()
	for (collection in collections) {
		collectionFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "collection.txt")
		if (file.exists(collectionFile)) {
			message("\tIn '", collection, "', found collection annotation file:", collectionFile);
			collectionDT = fread(collectionFile);
			setnames(collectionDT, tolower(colnames(collectionDT)));
			missCols = setdiff(collectionColNames, colnames(collectionDT));
			for (col in missCols) collectionDT[, col:=NA, with=FALSE];
			collectionDT = collectionDT[,collectionColNames, with=FALSE] #subset

		} else {
			message("\tIn '", collection, "', no collection file.");
			collectionDT = as.data.table(setNames(replicate(length(collectionColNames), NA, simplify = FALSE), collectionColNames));
		}
		collectionDT[,collectionname:=collection]
		collectionsDT = rbind(collectionsDT, collectionDT);
	}
	setkey(collectionsDT, "collectionname")
	setcolorder(collectionsDT, c("collectionname", collectionColNames));
	collectionsDT
}

#
#' Given a folder containing region collections in subfolders, this function
#' will either read the annotation file if one exists, or create a generic
#' annotation file.

#' @param dbLocation	folder where your regionDB is stored.
#' @param filePattern	passed to list.files; you can use this
#'	to select only certain file names in your folders.
#' @param refreshSizes	should I recreate the sizes files 
#'	documenting how many regions (lines) are in each region set?
#' 
#' @return Region set annotation (data.table)
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionAnno = readRegionSetAnnotation(dbLocation=dbPath)
readRegionSetAnnotation = function(dbLocation, 
					filePattern = "", 
					refreshSizes=FALSE) {
	size=NULL # Silence R CMD check Notes.
	#Build a data.table annotating the beds.
	#Should give collections
	collections = list.dirs(path=dbLocation, full.names=FALSE, recursive=FALSE)
	message("Reading region annotations...");
	annoDT = data.table();
	# Define pre-approved column names (others will be ignored)
	annotationColNames = c("filename", "cellType", "description", "tissue", "dataSource", "antibody", "treatment")

	for (collection in collections) {
		files = list.files(paste0(dbLocation,"/",collection, "/regions"), filePattern)
		#eliminate special annotation files
		#specialFileInd = grep("^0", files)
		#if (length(specialFileInd) > 0) {
		#	files = files [ -specialFileInd]
		#}
		if (length(files) <1) { 
			message("\tIn '", collection, "', no files found.");
			next;
		}
		collectionAnnoDT = data.table(collection=collection, filename=files, size=-1); #preserve new ones
		setkey(collectionAnnoDT, "filename")

		#look for index file
		indexFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "index.txt")
		if (file.exists(indexFile)) {
			message("\tIn '", collection, "', found index file:", indexFile);
			indexDT = fread(indexFile);
			indexDT[,filename:=as.character(filename)]
			setnames(indexDT, tolower(colnames(indexDT)));
		} else {
			message("\tIn '", collection, "', no index file. Found ", length(files), " files to load with defaults (filename only)");
#			indexDT = as.data.table(setNames(replicate(length(annotationColNames),character(0), simplify = F), annotationColNames));
			indexDT= data.table(filename=files)
		}
			missCols = setdiff(tolower(annotationColNames), colnames(indexDT));
			for (col in missCols) indexDT[, col:=NA, with=FALSE];
			indexDT = indexDT[,tolower(annotationColNames), with=FALSE]
 #subset
			# Revert back to camelCase
			setnames(indexDT, tolower(annotationColNames), annotationColNames);
		setkey(indexDT, "filename")
		#look for size file
		sizeFile = paste0(enforceTrailingSlash(dbLocation), enforceTrailingSlash(collection), "sizes.txt")
		if (file.exists(sizeFile) & !refreshSizes) {
			groupSizes = fread(sizeFile)
			#collectionAnnoDT[,size:=-1]
			setkey(groupSizes, "filename");
			groupSizes[, size_int:=as.double(size)]
			collectionAnnoDT[groupSizes, size:=size_int]
		}
		if (any(collectionAnnoDT[,size] < 0)) {
			message("Collection: ", collection, ". Creating size file...")
			collectionAnnoDT[,size:=countFileLines(paste0(dbLocation, "/", collection, "/regions/", filename)), by=filename]
			write.table(collectionAnnoDT[,list(filename,size)], file=sizeFile, quote=FALSE, row.names=FALSE, sep="\t")
		}

	collectionAnnoDT = indexDT[collectionAnnoDT]
#	collectionAnnoDT = collectionAnnoDT[indexDT]
	annoDT = rbind(annoDT, collectionAnnoDT);
	} #end loop through collections

	return(annoDT)
}



#' This function takes a region annotation object and reads in the regions,
#' returning a GRangesList object of the regions.
#' 
#' @param dbLocation	folder of regiondB
#' @param annoDT	output of readRegionSetAnnotation().
#' @param useCache	uses simpleCache to cache and load the results
#' @param limit	for testing purposes, limit the nmber of files read.
#'	NULL for no limit (default).
#'
#' @return GRangesList object
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionAnno = readRegionSetAnnotation(dbLocation=dbPath)
#' regionGRL = readRegionGRL(dbLocation= dbPath, regionAnno, useCache=FALSE)
readRegionGRL = function(dbLocation, annoDT, useCache=TRUE, limit=NULL) {
	grl = GRangesList()
	dbLocation = enforceTrailingSlash(dbLocation);
	
	for (iCol in unique(annoDT$collection)) {
	message(iCol);
	filesToRead = annoDT[collection==iCol,list(fullFilename=paste0(dbLocation, sapply(collection, enforceTrailingSlash), "regions/", filename)), by=filename]$fullFilename
	if (useCache) {
		if (requireNamespace("simpleCache", quietly=TRUE)) {
			simpleCache::simpleCache(iCol, {readCollection(filesToRead)}, cacheDir=paste0(dbLocation, iCol), buildEnvir=nlist(filesToRead))
		} else {
			warning("You don't have simpleCache installed, so you won't be able to cache the regionDB after reading it in. Install simpleCache to speed up later database loading.")
		}
	} else {
	assign(iCol, readCollection(filesToRead, limit));
	}
	grl = c(grl, get(iCol))
	}
	return(grl)
}

#' Given a bunch of region set files, read in all those flat (bed) files and create a
#' GRangesList object holding all the region sets. This function is used by readRegionGRL
#' to process annotation objects.
#'
#' @param filesToRead	a vector containing bed files
#' @param limit	for testing purposes, limit the number of files read. NULL for no limit (default).
#'
#' @export
#' @examples
#' files = list.files(system.file("extdata", "hg19/ucsc_example/regions",
#'	 package="LOLA"), pattern="*.bed")
#' regionAnno = readCollection(files)
readCollection = function(filesToRead, limit=NULL) {
	grl = GRangesList()
	if (is.null(limit)) {
		limit = length(filesToRead);
	}	else {
		message("limit files: ", limit);
	}
	message("Reading ", length(filesToRead), " files...")
	for (i in 1:limit) {
		message(i, ": ", filesToRead[i]);
		filename = filesToRead[i]
		if (file.exists(filename)) {
			success = tryCatch( { 
				DT = fread(paste0(filename))
				cn = colnames(DT)
				cn[1]
				tfbsgr = dtToGr(DT, cn[1], cn[2], cn[3]);
				grl[[i]] = tfbsgr;
				TRUE
			},
			error = function(e) { message(i, " ERR:", filename); print(e); return(FALSE); } )
			if (!success) { grl[[i]] = GRanges(); }
		} else {
			message("Skipping (file not found):", filename);
			grl[[i]] = GRanges();
		}
	}
	return(grl);
}


#' Given two regionDBs, (lists returned from readRegionDB()),
#' This function will combine them into a single regionDB. This
#' will enable you to combine, for example, LOLA Core databases
#' with custom databases into a single analysis.
#'
#' @param dbA First regionDB database.
#' @param dbB Second regionDB database.
#'
#' @return A combined regionDB.
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionDB = loadRegionDB(dbPath)
#' combinedRegionDB = mergeRegionDBs(regionDB, regionDB)

mergeRegionDBs = function(dbA, dbB) {
	myNames = names(regionDB)
	# Loop through each item and concat them
	combinedRegionDB = list()	
	for (item in myNames) {
		message(item)
		if ("character" %in% class(dbA[[item]])) {
			combinedRegionDB[[item]] = c(dbA[[item]], dbB[[item]])
		} else if ("data.table" %in% class(dbA[[item]])) {
			combinedRegionDB[[item]] = rbind(dbA[[item]], dbB[[item]])
		} else if ("GRangesList" %in% class(dbA[[item]])) {
			combinedRegionDB[[item]] = c(dbA[[item]], dbB[[item]])
		}			
	}
	return(combinedRegionDB)
}
