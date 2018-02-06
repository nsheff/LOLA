######################################################################
# Generic Region DB Loading Functions

######################################################################


#' Helper function to annotate and load a regionDB, a folder with
#' subfolder collections of regions.
#'
#' @param dbLocation	folder where your regionDB is stored, or list of such folders
#' @param useCache	uses simpleCache to cache and load the results
#' @param limit 	You can limit the number of regions for testing.
#'	Default: NULL (no limit)
#' @param collections Restrict the database loading to this list of collections
#'
#' @return regionDB list containing database location, region and
#' collection annotations, and regions GRangesList
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionDB = loadRegionDB(dbLocation=dbPath)

loadRegionDB = function(dbLocation, useCache=TRUE, limit=NULL,
collections=NULL) {
	# Base case
	if (length(dbLocation) == 1) {
		collectionAnno = readCollectionAnnotation(dbLocation, collections)
		regionAnno = readRegionSetAnnotation(dbLocation, collections)
		regionGRL = readRegionGRL(dbLocation, regionAnno, useCache=useCache, limit=limit)
		return(nlist(dbLocation, regionAnno, collectionAnno, regionGRL))
	} else {
		return(
			mergeRegionDBs(
				loadRegionDB(dbLocation[1]),
				loadRegionDB(dbLocation[-1])
			)
		)
	}
}

#' Read collection annotation
#'
#' @param dbLocation	Location of the database
#' @param collections Restrict the database loading to this list of collections.
#' Leave NULL to load the entire database (Default).
#'
#' @return Collection annotation data.table
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' collectionAnno = readCollectionAnnotation(dbLocation=dbPath)
readCollectionAnnotation = function(dbLocation, collections=NULL) {
	annoDT = data.table()
	collectionList = list.dirs(path=dbLocation, full.names=FALSE, recursive=FALSE)
	if (! is.null(collections)) {
		# Restrict the list if parameter is passed.
		collectionList = intersect(collections, collectionList)
	}
	if (length(collectionList) == 0) {
		stop(paste0("No collections in ", dbLocation, ". Check your path."))
	}
	message("Reading collection annotations: ", paste(collections, collapse=", "))
	collectionColNames = c("collector", "date", "source", "description")
	collectionsDT = data.table()
	for (collection in collectionList) {
		collectionFile = paste0(enforceTrailingSlash(dbLocation),
enforceTrailingSlash(collection), "collection.txt")
		if (file.exists(collectionFile)) {
			message("\t", collection, ": found collection annotation:", collectionFile)
			collectionDT = fread(collectionFile)
			setnames(collectionDT, tolower(colnames(collectionDT)))
			missCols = setdiff(collectionColNames, colnames(collectionDT))
			for (col in missCols) collectionDT[, (col) := NA]
			collectionDT = collectionDT[, collectionColNames, with=FALSE] #subset
		} else {
			regionFiles = list.files(paste0(dbLocation,"/",collection, "/regions"))
			if (length(regionFiles) < 1) {
				warning(collection, " has no collection annotation or region files. Skipping...")
				next
			} else {
			message("\tIn collection '", collection,
			 "', consider adding a 'collection.txt' annotation file.")
			collectionDT = as.data.table(
				setNames(replicate(length(collectionColNames), NA, simplify = FALSE),
				collectionColNames))
			}
		}

		collectionDT[,collectionname:=collection]
		collectionsDT = rbind(collectionsDT, collectionDT)
	}
	if (nrow(collectionsDT) < 1) {
		stop("No regions found. Are you sure '", dbLocation, "' is a region database?")
	}
	setkey(collectionsDT, "collectionname")
	setcolorder(collectionsDT, c("collectionname", collectionColNames))
	collectionsDT
}

#
#' Given a folder containing region collections in subfolders, this function
#' will either read the annotation file if one exists, or create a generic
#' annotation file.

#' @param dbLocation	folder where your regionDB is stored.
#' @param collections Restrict the database loading to this list of collections
#' Leave NULL to load the entire database (Default).
#' @param refreshCaches	should I recreate the caches? Default: FALSE
#' @param refreshSizes should I refresh the size files? Default:TRUE
#' @param useCache Use simpleCache to store results and load them?
#'
#' @return Region set annotation (data.table)
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionAnno = readRegionSetAnnotation(dbLocation=dbPath)
readRegionSetAnnotation = function(dbLocation, collections = NULL,
					refreshCaches=FALSE, refreshSizes=TRUE,
					useCache=TRUE) {
	size=NULL # Silence R CMD check Notes.
	# Build a data.table annotating the beds.
	# Should give collections
	annoDT = data.table()
	dbLocation = enforceTrailingSlash(dbLocation)
	collectionList = list.dirs(path=dbLocation, full.names=FALSE, recursive=FALSE)
	if (! is.null(collections)) {
		# Restrict the list if parameter is passed.
		collectionList = intersect(collections, collectionList)
	}
	message("Reading region annotations...")
	if (length(collectionList) == 0) {
		stop(paste0("No collections were found in ", dbLocation, ". Check your path."))
	}
	for (collection in collectionList) {
		if (useCache & requireNamespace("simpleCache", quietly=TRUE)) {
			simpleCache::simpleCache(paste0(collection, "_files"), {
				readCollectionFiles(dbLocation, collection, refreshSizes=refreshSizes)},
				cacheDir=enforceTrailingSlash(paste0(dbLocation, collection)),
				buildEnvir =nlist(dbLocation, collection), recreate=refreshCaches)
		} else {
			warning("You don't have simpleCache installed, so you won't be able to cache the
regionDB after reading it in. Install simpleCache to speed up later database
loading.")
			assign(paste0(collection, "_files"),
			readCollectionFiles(dbLocation, collection, refreshSizes=refreshSizes))
		}

		annoDT = rbind(annoDT, get(paste0(collection, "_files")))
	} #end loop through collections


	if (nrow(annoDT) < 1) {
		stop("No regions found. Are you sure '", dbLocation, "' is a region database?")
	}
	return(annoDT)
}
#' Given a database and a collection, this will create the region annotation
#' data.table; either giving a generic table based on file names, or by
#' reading in the annotation data.
#' @param dbLocation	folder where your regionDB is stored.
#' @param collection Collection folder to load
#' @param refreshSizes	should I recreate the sizes files
#'	documenting how many regions (lines) are in each region set?
#' @return A data.table annotating the regions in the collections.
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionAnno = readCollectionFiles(dbLocation=dbPath, "ucsc_example")
readCollectionFiles = function(dbLocation, collection, refreshSizes=FALSE) {
	# Define pre-approved column names (others will be ignored)
	annotationColNames = c("filename", "cellType", "description", "tissue",
	"dataSource", "antibody", "treatment")
	message(paste0(dbLocation,"/",collection, "/regions"))
	files = list.files(paste0(dbLocation,"/",collection, "/regions"))
		if (length(files) <1) {
			message("\tIn '", collection, "', no files found.")
			return()
		}

		# Preserve new ones
		collectionAnnoDT = data.table(collection=collection, filename=files, size=-1)
		setkey(collectionAnnoDT, "filename")

		# Look for index file
		indexFile = paste0(
			enforceTrailingSlash(dbLocation),
			enforceTrailingSlash(collection),
			"index.txt")
		if (file.exists(indexFile)) {
			message("\tIn '", collection, "', found index file:", indexFile)
			indexDT = fread(indexFile)
			indexDT[, filename:=as.character(filename)]
			setnames(indexDT, tolower(colnames(indexDT)))
		} else {
			message("\tIn '", collection, "', no index file. Found ", length(files), " files
to load with defaults (filename only)")
			indexDT = data.table(filename=files)
		}
			missCols = setdiff(tolower(annotationColNames), colnames(indexDT))

			# Populate any missing columns with NAs (of character type):
			for (col in missCols) indexDT[, (col) := as.character(NA)]
			indexDT = indexDT[,tolower(annotationColNames), with=FALSE]
			# Subset
			# Revert back to camelCase
			setnames(indexDT, tolower(annotationColNames), annotationColNames)
		setkey(indexDT, "filename")
		#look for size file
		sizeFile = paste0(
			enforceTrailingSlash(dbLocation),
			enforceTrailingSlash(collection),
			"sizes.txt")
		if (file.exists(sizeFile) & !refreshSizes) {
			groupSizes = fread(sizeFile)
			#collectionAnnoDT[,size:=-1]
			setkey(groupSizes, "filename")
			groupSizes[, size_int:=as.double(size)]
			collectionAnnoDT[groupSizes, size:=size_int]
		}
		if (any(collectionAnnoDT[,size] < 0)) {
			message("Collection: ", collection, ". Creating size file...")
			collectionAnnoDT[,
				size:=countFileLines(paste0(dbLocation, "/", collection, "/regions/", filename)),
				by=filename]
			write.table(collectionAnnoDT[,list(filename,size)],
				file=sizeFile, quote=FALSE,
				row.names=FALSE, sep="\t")
		}

	collectionAnnoDT = indexDT[collectionAnnoDT]

	# Let's try to avoid NAs in the description columns:
	buildGenericDescription = function(cellType, tissue, antibody, treatment, collection) {
		description = collection
		# First, try either a generic cellType or tissue
		if (any(!is.na(cellType))) {
			description[!is.na(cellType)] =
				paste(description[!is.na(cellType)], cellType[!is.na(cellType)])
		} else if (any(!is.na(tissue))) {
			description[!is.na(tissue)] =
				paste(description[!is.na(tissue)], tissue[!is.na(tissue)])
		}

		# Second, try either a generic antibody or treatment
		if (any(!is.na(antibody))) {
			description[!is.na(antibody)] =
				paste(description[!is.na(antibody)], antibody[!is.na(antibody)])
		} else if (any(!is.na(treatment))) {
			description[!is.na(treatment)] =
				paste(description[!is.na(treatment)], treatment[!is.na(treatment)])
		}
		return(description)
	}
	collectionAnnoDT[is.na(description),
description:=buildGenericDescription(cellType, tissue, antibody, treatment,
collection)]

	return(collectionAnnoDT)
}


#' This function takes a region annotation object and reads in the regions,
#' returning a GRangesList object of the regions.
#'
#' @param dbLocation	folder of regiondB
#' @param annoDT	output of readRegionSetAnnotation().
#' @param refreshCaches	should I recreate the caches?
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
readRegionGRL = function(dbLocation, annoDT, refreshCaches=FALSE,
	useCache=TRUE, limit=NULL) {
	grl = GRangesList()
	dbLocation = enforceTrailingSlash(dbLocation)

	for (iCol in unique(annoDT$collection)) {
	message(iCol)
	filesToRead = annoDT[collection==iCol,list(fullFilename=paste0(dbLocation,
sapply(collection, enforceTrailingSlash), "regions/", filename)),
by=filename]$fullFilename
	if (useCache & requireNamespace("simpleCache", quietly=TRUE)) {
		simpleCache::simpleCache(iCol, {
			readCollection(filesToRead)},
			cacheDir=paste0(dbLocation, iCol),
			buildEnvir=nlist(filesToRead), recreate=refreshCaches)
	} else {
		if (length(filesToRead) > 100) {
			# Notify user of caching possibility on large databases
			message("Install simpleCache to speed up future database loading.")
		}
		assign(iCol, readCollection(filesToRead, limit))
	}
	grl = c(grl, get(iCol))
	}
	return(grl)
}

#' Given a bunch of region set files, read in all those flat (bed) files and
#' create a GRangesList object holding all the region sets. This function is
#' used by readRegionGRL to process annotation objects.
#'
#' @param filesToRead	a vector containing bed files
#' @param limit	for testing purposes, limit the number of files read.
#' NULL for no limit (default).
#'
#' @return A GRangesList with the GRanges in the filesToRead.
#' @export
#' @examples
#' files = list.files(system.file("extdata", "hg19/ucsc_example/regions",
#'	 package="LOLA"), pattern="*.bed")
#' regionAnno = readCollection(files)
readCollection = function(filesToRead, limit=NULL) {
	grl = GRangesList()
	if (is.null(limit)) {
		limit = length(filesToRead)
	}	else {
		message("limit files: ", limit)
	}
	message("Reading ", length(filesToRead), " files...")
	for (i in seq_len(limit)) {
		message(i, ": ", filesToRead[i])
		filename = filesToRead[i]
		if (file.exists(filename)) {
			success = tryCatch( {
				tfbsgr = readBed(filename)
				# DT = fread(paste0(filename))
				# cn = colnames(DT)
				# cn[1]
				# tfbsgr = dtToGr(DT, cn[1], cn[2], cn[3])
				grl[[i]] = tfbsgr
				TRUE
			},
			error = function(e) {
				message(i, " ERR:", filename)
				message(e)
				return(FALSE)
			})
			if (!success) { grl[[i]] = GRanges() }
		} else {
			message("Skipping (file not found):", filename)
			grl[[i]] = GRanges()
		}
	}
	return(grl)
}


#' Given two regionDBs, (lists returned from loadRegionDB()),
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
	myNames = names(dbA)
	# Loop through each item and concat them
	combinedRegionDB = list()
	for (item in myNames) {
		message(item)
		if ("character" %in% class(dbA[[item]])) {
			combinedRegionDB[[item]] = c(dbA[[item]], dbB[[item]])
		} else if ("data.table" %in% class(dbA[[item]])) {
			combinedRegionDB[[item]] = rbind(dbA[[item]], dbB[[item]])
		} else if (is(dbA[[item]], "GRangesList")) {
			combinedRegionDB[[item]] = c(dbA[[item]], dbB[[item]])
		}
	}
	return(combinedRegionDB)
}



#' Grab a single region set from a database, specified by filename.
#'
#' If you want to work with a LOLA regionDB region set individually,
#' this function can help you. It can extract individual (or subsets of)
#' region sets from either loaded regionDBs, loaded with loadRegionDB(), or
#' from a database on disk, where only the region sets of interest will
#' be loaded.
#'
#' @param regionDB A region database loaded with loadRegionDB().
#' @param filenames Filename(s) of a particular region set to grab.
#' @param collections (optional) subset of collections to list
#'
#' @return A GRanges object derived from the specified file in the regionDB.
#' @export
#' @example
#' R/examples/example.R
getRegionSet = function(regionDB, filenames, collections=NULL) {
	if ("regionAnno" %in% names(regionDB)) {
		# it's a loaded regionDB object, we just extract the region set.
	
		ind = regionDB$regionAnno[,which(filename %in% filenames)]
		if (length(ind) < 1) {
			stop("That filename was not found in the database.")
		}
		return(regionDB$regionGRL[ind])
	} else {
		filesToRead = getRegionFile(regionDB, filenames, collections)
		grl = readCollection(filesToRead)
		return(grl)
	}
}

#' Grab the filename for a a single region set from a database specified by filename.
#'
#' Like getRegionSet but returns a filename instead of a GRanges object. Given a local
#' filename, returns a complete absolulte path so you can read that file in.
#'
#' @param dbLocation	folder of regionDB
#' @param filenames Filename(s) of a particular region set to grab.
#' @param collections (optional) subset of collections to list
#'
#' @return A filename the specified file in the regionDB.
#' @export
#' @example
#' R/examples/example.R
getRegionFile = function(dbLocation, filenames, collections = NULL) {
	dbLocation = enforceTrailingSlash(dbLocation)
	collectionAnno =
		readCollectionAnnotation(dbLocation, collections)
	regionAnno =
		readRegionSetAnnotation(dbLocation, collections)
	filesToRead =
		regionAnno[filename %in% filenames, list(fullFilename = paste0(dbLocation,
		sapply(collection, enforceTrailingSlash), "regions/",
		filename)), by = filename]$fullFilename
	return(filesToRead)
}


#' Lists the region sets for given collection(s) in a region database on disk.
#'
#' @param regionDB File path to region database
#' @param collections (optional) subset of collections to list
#'
#' @return a list of files in the given collections
#' @export
#' @examples
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' listRegionSets(dbPath)
listRegionSets = function(regionDB, collections=NULL) {
	if (is.character(regionDB)){
		dbLocation = regionDB
		dbLocation = enforceTrailingSlash(dbLocation)
		collectionAnno = readCollectionAnnotation(dbLocation,
			collections)
		regionAnno = readRegionSetAnnotation(dbLocation, collections)
	} else if ("dbLocation" %in% names(regionDB)) {
		regionAnno = regionDB$regionAnno
	} else {
		stop("This is not a proper regionDB or path to one.")
	}
	return(regionAnno$filename)
}



