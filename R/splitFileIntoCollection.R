#filename = "~/fhgfs/share/regionDB/hg19/sheffield_dnase/TableS03-dhs-to-cluster.txt"
#splitCol="refined_cluster"
#' This function will take a single large bed file that is annotated
#' with a column grouping different sets of similar regions, 
#' and split it into separate files for use with the LOLA collection format.
#'
#' @param filename	the file to split
#' @param splitCol	factor column that groups the lines in the file by set
#'
#' @export
#' @examples
#' combFile = system.file("extdata", "examples/combined_regions.bed", package="LOLA")
#' splitFileIntoCollection(combFile, 4)
splitFileIntoCollection = function(filename, splitCol) {
	DT = fread(paste0(filename))
	collectionFolder = paste0(dirname(filename), "/", basename(filename), "_collection")
	dir.create(collectionFolder);
	sDT = splitDataTable(DT, splitFactor=splitCol)
	nDT = names(sDT);
	for (i in 1:length(sDT)) {
		message(nDT[[i]]);
		print(sDT[[i]]);
		write.tsv(sDT[[i]], paste0(collectionFolder, "/", nDT[[i]], ".bed"))
	}
	message("Collection written to ", collectionFolder)
}


