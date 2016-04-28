#filename = "~/fhgfs/share/regionDB/hg19/sheffield_dnase/TableS03-dhs-to-cluster.txt"
#splitCol="refined_cluster"
#' This function will take a single large bed file that is annotated
#' with a column grouping different sets of similar regions,
#' and split it into separate files for use with the LOLA collection format.
#'
#' @param filename	the file to split
#' @param splitCol	factor column that groups the lines in the file by set. It 
#' 	should be an integer.
#' @param collectionFolder	name of folder to place the new split files.
#' @param filenamePrepend	string to prepend to the filenames. Defaults to blank.
#'
#' @return No return value.
#' @export
#' @examples
#' combFile = system.file("extdata", "examples/combined_regions.bed", package="LOLA")
#' splitFileIntoCollection(combFile, 4)
splitFileIntoCollection = function(filename, splitCol, collectionFolder=NULL,
		filenamePrepend="") {
	DT = fread(paste0(filename))
	if (is.null(collectionFolder)) {
		# Default collection folder:
		collectionFolder = paste0(dirname(filename), "/", basename(filename), "_collection")
	}
	dir.create(collectionFolder)
	sDT = splitDataTable(DT, splitFactor=splitCol)
	nDT = names(sDT)
	nDT = paste0(filenamePrepend, nDT)
	for (i in seq_along(sDT)) {
		message(nDT[[i]])
		write.tsv(sDT[[i]], paste0(collectionFolder, "/", nDT[[i]], ".bed"), col.names=FALSE)
	}
	message("Collection written to ", collectionFolder)
}
