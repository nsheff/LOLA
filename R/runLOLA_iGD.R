######################################################################
# ENRICHMENT - Actual workhorse enrichment calculation functions
######################################################################

#' Enrichment Calculation
#'
#' Workhorse function that uses iGD to calculate overlaps between userSets,
#' and then uses a fisher's exact test rank them by significance
#' of the overlap.
#'
#' @param userSets		Regions of interest (Can be GRanges objects or a list)
#' @param userUniverse	Regions tested for inclusion in userSets (Can be a GRanges object)
#' @param pepRegionDB	Region DB to check for overlap, output list object from loadPEPdb()
#' @param cores	Number of processors
#' @param redefineUserSets	run redefineUserSets() on your userSets?
#' @param direction    Defaults to "enrichment", but may also accept
#'     "depletion", which will swap the direction of the fisher test (use
#'     'greater' or less' value passed to the 'alternative' option of
#'     fisher.test)
#' @return Data.table with enrichment results. Rows correspond to individual
#' pairwise fisher's tests comparing a single userSet with a single databaseSet.
#' The columns in this data.table are: userSet and dbSet: index into their
#' respective input region sets. pvalueLog: -log10(pvalue) from the fisher's exact
#' result; oddsRatio: result from the fisher's exact test; support: number of
#' regions in userSet overlapping databaseSet; rnkPV, rnkOR, rnkSup: rank in this
#' table of p-value, oddsRatio, and Support respectively. The --value is the
#' negative natural log of the p-value returned from a one-sided fisher's exact
#' test. maxRnk, meanRnk: max and mean of the 3 previous ranks, providing a
#' combined ranking system. b, c, d: 3 other values completing the 2x2 contingency
#' table (with support). The remaining columns describe the dbSet for the row.
#'
#' If you have the qvalue package installed from bioconductor, runLOLA will add
#' a q-value transformation to provide FDR scores automatically.
#' @export
#' @example
#' R/examples/example.R


setA = LOLA::readBed("~/lola_vignette_data/setA_complete.bed")
setB = LOLA::readBed("~/lola_vignette_data/setB_complete.bed")
setC = LOLA::readBed("~/lola_vignette_data/setC_complete.bed")
active_dhs = LOLA::readBed("~/lola_vignette_data/activeDHS_universe.bed")
usersets = GRangesList(setA, setB, setC)
pepregiondb = loadPEPdb("/project/shefflab/resources/regions/LOLAHema/hg38/test_bedset4/test_bedset4_PEP/test_bedset4_config.yaml")

# Set igd db for testing purposes
igd_obj = IGDr::IGDr(pepregiondb[[5]])
setA_df = as.data.frame(setA)
ch_names = setA_df$seqnames
start = setA_df$start
end = setA_df$end

# Set universe vectors for testing
universe_df = as.data.frame(active_dhs, row.names=seq_along(active_dhs))
univ_names = universe_df$seqnames
univ_start = universe_df$start
univ_end = universe_df$end
igd_univ = IGDr::search_nr(igd_obj, nrow(universe_df), univ_names, univ_start, univ_end)

# Compare igd search_nr vs countOverlaps
igd_ol = IGDr::search_nr(igd_obj, length(setA), ch_names, start, end)
gr_ol = countOverlaps(pepregiondb[[4]], setA, type = c("any"), minoverlap = 0)
gr_univ = countOverlaps(pepregiondb[[4]], active_dhs, )

# Test runLOLA with iGD vs original function using countOverlaps
newLOLA_res =  runLOLA2(usersets, active_dhs, pepregiondb, cores = 1, direction = "enrichment") # time elapsed=1.959 sec
oldLOLA_res = runLOLA(usersets, active_dhs, pepregiondb, cores = 1, direction = "enrichment") # time elapsed=4.885 sec


####### Define new runLOLA function
runLOLA2 = function(userSets, userUniverse, pepRegionDB, cores=1,
                   redefineUserSets=FALSE, direction="enrichment") {
  
  if(!(is(userSets, "GRanges") || is(userSets, "GRangesList"))) {
    stop("userSets should be a GRanges object or GRanges list. Check object class")
  }
  if(!(is(userUniverse, "GRanges")) {
    stop("userUniverse should be a GRanges object. Check object class")
  }
  
  
  # Silence R CMD check Notes:
  support=d=b=userSet=pValueLog=rnkSup=rnkPV=rnkOR=NULL
  oddsRatio=maxRnk=meanRnk=dbSet=description=NULL
  annotationDT = pepRegionDB$regionAnno
  IGDreferenceLoc = pepRegionDB$iGDRefDatabase
  IGDindex = pepRegionDB$iGDRefIndex
  
  if (direction == "depletion") {
    fisherAlternative = "less"
  } else {
    fisherAlternative = "greater"
  }
  
  annotationDT[, dbSet := seq_len(nrow(annotationDT))]
  setkey(annotationDT, dbSet)
  
  ### Data sanity checks ###
  #Confirm we received GRangesList objects, convert from list if possible.
  userSets = GRangesList(userSets)
  #setLapplyAlias(cores)
  
  
  #if (any(is.null(names(testSetsGRL)))) {
    #names(testSetsGRL) = seq_along(testSetsGRL)
  #}
  
  if (redefineUserSets) { #redefine user sets in terms of universe?
    userSets =	redefineUserSets(userSets, userUniverse, cores=cores)
    userSets = listToGRangesList(userSets)
  }
  
  if (! any( isDisjoint( userSets) ) ) {
    message("You have non-disjoint userSets.")
  }
  
  ### Construct significance tests ###
  message("Calculating unit set overlaps...")
  
  
  # Returns for each userSet, a vector of length length(testSetsGRL), with total
  # number of regions in that set overlapping anything in each testSetsGRL; this
  # is then lapplied across each userSet.
  # ----------------------
  # Replace countoverlaps function with igd search
  # Need to convert refdatabase.igd into an IGDr object(open/load an igd database for search )
  
  #geneSetDatabaseOverlap =
  #lapplyAlias( (userSets), countOverlapsRev, testSetsGRL, minoverlap=minOverlap)
  
  
  # Convert usersets GRanges obj to a data frame for overlaps calculation
  
  #userSetsData = as.data.frame(userSets)
  #userSetsLength = nrow(userSetsData)
  #chromNames = userSetsData$seqnames
  #chromStart = userSetsData$start
  #chromEnd = userSetsData$end
  # Get the number of queries to be searched
  #queriesN = nrow(userSetsData)
  
  # Produce a vector to match order of overlaps with order of files in annotation
  igdIndexFiles = IGDindex$File
  annotation = pepRegionDB$regionAnno
  AnnoFileNames = sapply(annotation$output_file_path, basename)
  correctRefOrder = match(AnnoFileNames, igdIndexFiles)
  
  # Alternative approach  
  userSetsLength = unlist(lapply((userSets), length))
  userSetsData = lapply(userSets, as.data.frame)
  IGDrefDB = IGDr::IGDr(IGDreferenceLoc)
  IGDoverlapList = list()
  for (i in seq_along(userSetsData)) {
    chroms = userSetsData[[i]]$seqnames
    starts = userSetsData[[i]]$start
    ends = userSetsData[[i]]$end
    userSetLength = nrow(userSetsData[[i]])
    IGDoverlapList[[i]] = IGDr::search_nr(IGDrefDB, userSetLength, chroms, starts, ends)
    IGDoverlapList[[i]] = IGDoverlapList[[i]][correctRefOrder]
  }
  

  # Convert iGD db to an IGDr object and perform the overlap calculation
  #IGDrefDB = IGDr::IGDr(IGDreferenceLoc)
  
  #IGDrefDatabaseOverlap = IGDr::search_nr(IGDrefDB, queriesN, 
                                          #chromNames, chromStart, chromEnd)

  # This will become "support" -- the number of regions in the
  # userSet (which I implicitly assume is ALSO the number of regions
  # in the universe) that overlap anything in each database set.
  # Turn results into an overlap matrix. It is
  # dbSets (rows) by userSets (columns), counting overlap.
  
  #overlapMatrix = do.call(cbind, IGDrefDatabaseOverlap)
  IGDoverlapMatrix = do.call(cbind, IGDoverlapList)
  
  message("Calculating universe set overlaps...")
  # Now for each test set, how many items *in the universe* does
  # it overlap? This will go into the calculation for c
  
  #faster. Returns number of items in userUniverse.
  
  #testSetsOverlapUniverse = countOverlaps(testSetsGRL, userUniverse, minoverlap=minOverlap)
  
  # Convert useruniverse into a df and perform overlap calculations with igd ref database
  
  universeData = as.data.frame(userUniverse, row.names=seq_along(userUniverse)) # universe names req below
  universeNames = universeData$seqnames
  universeStart = universeData$start
  universeEnd = universeData$end
  universeRegionN = nrow(universeData)
  
  IGDuniverseOverlap = IGDr::search_nr(IGDrefDB, universeRegionN, 
                                       universeNames, universeStart, universeEnd)
  IGDuniverseOverlap = IGDuniverseOverlap[correctRefOrder]
  names(IGDuniverseOverlap) = seq_along(IGDuniverseOverlap)
  
  # Returns number of items in test set (not used:)
  #testSetsOverlapUniverse = countOverlapsAny(testSetsGRL, userUniverse)
  # Total size of the universe (could get size of either original GRanges or df since only 1 file)
  universeLength = length(userUniverse)
  
  # To build the fisher matrix, support is 'a'
  
  scoreTable = data.table(reshape2::melt(t(IGDoverlapMatrix), variable.factor=FALSE))

  setnames(scoreTable, c("Var1", "Var2", "value"), c("userSet", "dbSet", "support"))
  
  # reshape2 has an annoying habit of converting strings into factors, which
  # is undesirable. If the userSets are named with strings, make sure they stay
  # character. Integers are already handled appropriately.
  
  if ("factor" %in% class(scoreTable[, userSet])) {
    scoreTable$userSet = as.character(scoreTable$userSet)
  }
  
  
  message("Calculating Fisher scores...")
  # b = the # of items *in the universe* (1 bed file) that overlap each dbSet (bed file),
  # less the support; This is the number of items in the universe
  # that are in the dbSet ONLY (not in userSet)
  # c = the size of userSet, less the support; This is the opposite:
  # Items in the userSet ONLY (not in the dbSet)
  
  scoreTable[,c("b", "c"):=list(
    b=IGDuniverseOverlap[match(dbSet, names(IGDuniverseOverlap))]-support, c=userSetsLength-support)]
  
  
  # This is the regions in the universe, but not in dbSet nor userSet.
  scoreTable[,d:=universeLength-support-b-c]
  if( scoreTable[,any(b<0)] ) { # Inappropriate universe.
    warning(cleanws("Negative b entry in table. This means either: 1) Your user sets
                    contain items outside your universe; or 2) your universe has a region that
                    overlaps multiple user set regions, interfering with the universe set overlap
                    calculation."))
    
    return(scoreTable)
  }
  if( scoreTable[,any(c<0)] ) {
    warning("Negative c entry in table. Bug with userSetsLength; this should not happen.")
    return(scoreTable)
  }

  # Code for Fisher test looks fine as it is.
  scoreTable[,c("pValueLog", "oddsRatio") := 
               fisher.test(matrix(c(support,b,c,d), 2, 2), alternative=fisherAlternative)[c("p.value", "estimate")], by=list(userSet,dbSet)]
  
  # Include qvalue if package exists.
  if (requireNamespace("qvalue", quietly=TRUE)) {
    # Wrap in try block since this is not vital.
    # if you want qvalues...
    tryCatch( {
      scoreTable[,qValue:=qvalue::qvalue(pValueLog)$qvalue]
    }, error = function(e) { warning("Problem in FDR calculation with qvalue.") })
  } else {
    # Another possibility for the future:
    # scoreTable[,qValue:=qValues = pmin(pValues*length(pValues),1)]
  }
  scoreTable[, pValueLog:=-log10(pValueLog + 10^-322)]
  ### Finalize and Rank results ###
  scoreTable[, rnkSup:=rank(-support, ties.method="min"), by=userSet]
  scoreTable[, rnkPV:=rank(-pValueLog, ties.method="min"), by=userSet]
  scoreTable[, rnkOR:=rank(-oddsRatio, ties.method="min"), by=userSet]
  scoreTable[, maxRnk:=max(c(rnkSup, rnkPV, rnkOR)), by=list(userSet,dbSet)]
  scoreTable[, meanRnk:=signif(mean(c(rnkSup, rnkPV, rnkOR)), 3), by=list(userSet,dbSet)]
  
  # Append description column
  setkeyv(scoreTable, "dbSet")
  scoreTable = scoreTable[annotationDT]
  
  # limit description to 80 characters
  scoreTable[,description:=substr(description, 0, 80)]
  
  orderedCols = c("userSet", "dbSet", "collection", "pValueLog", "oddsRatio",
                  "support", "rnkPV", "rnkOR", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d",
                  "description", "cellType", "tissue", "antibody", "treatment", "dataSource", "filename")
  unorderedCols = setdiff(colnames(scoreTable), orderedCols)
  
  setcolorder(scoreTable, c(orderedCols, unorderedCols))
  
  scoreTable[order(pValueLog, -meanRnk, decreasing=TRUE),]
}









