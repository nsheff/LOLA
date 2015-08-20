#' A reduced GRanges object from the example regionDB database
#'
#'
#' @docType data
#' @keywords datasets
#' @name userUniverse
#' @usage data(sample_universe)
#' @format A GRanges object
#' @return No return value.
NULL

#' An example set of regions, sampled from the example database.
#'
#' A dataset containing a few sample regions.
#'
#' @docType data
#' @keywords datasets
#' @name userSet
#' @usage data(sample_input)
#' @format A GRanges object
#' @return No return value.
NULL

# This is how I produced the sample data sets:
#dbPath = system.file("extdata", "hg19", package="LOLA")
#regionDB = loadRegionDB(dbLocation= dbPath)
#userSet = reduce(do.call(c, (sampleGRL(regionDB$regionGRL, prop=c(.1,.25,.05,.05)))))
#userUniverse = reduce(do.call(c, regionDB$regionGRL))
#save(userSet, file="sample_input.RData")
#save(userUniverse, file="sample_universe.RData")
