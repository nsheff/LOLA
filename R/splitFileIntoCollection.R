
#filename = "~/fhgfs/share/regionDB/hg19/sheffield_dnase/TableS03-dhs-to-cluster.txt"
#splitCol="refined_cluster"
#' This function will take a single large bed file and split it into files
#' for use with the LOLA collection format.
#' @export
#'
splitFileIntoCollection = function(filename, splitCol) {
	DT = fread(paste0(filename))
	collectionFolder = paste0(dirname(filename), "/", basename(filename), "_collection")
	dir.create(collectionFolder);
	sDT = splitDataTable(DT, splitFactor=get("splitCol"))
	nDT = names(sDT);
	for (i in 1:length(sDT)) {
		message(nDT[[i]]);
		write.tsv(sDT[[i]], paste0(collectionFolder, "/", nDT[[i]], ".bed"))
	}		
}


