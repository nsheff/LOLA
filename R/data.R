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
#' @name userSets
#' @usage data(sample_input)
#' @format A GRangesList object
#' @return No return value.
#' @examples
#' \dontrun{
#'  This is how I produced the sample data sets:
#' dbPath = system.file("extdata", "hg19", package="LOLA")
#' regionDB = loadRegionDB(dbLocation= dbPath)
#' userSetA = reduce(do.call(c, (sampleGRL(regionDB$regionGRL,
#' prop=c(.1,.25,.05,.05,0)))))
#' userSetB = reduce(do.call(c, (sampleGRL(regionDB$regionGRL,
#' prop=c(.2,.05,.05,.05,0)))))
#'
#' userSets = GRangesList(setA=userSetA, setB=userSetB)
#' userUniverse = reduce(do.call(c, regionDB$regionGRL))
#' save(userSets, file="sample_input.RData")
#' save(userUniverse, file="sample_universe.RData")
#' }
NULL
