# LOLA: Genomic Locus Overlap Enrichment Analysis

The public-facing website for LOLA can be found at [http://databio.org/lola](http://databio.org/lola).

LOLA is an R package providing functions for testing overlap of sets of genomic regions with public and custom databases. You can think of it as testing your `bed file` (genome regions of interest) against a database of other `bed files` (regions from various previous studies) to look for enrichment of overlaps. This enables you to draw connections between newly generated data, and the growing public databases, leading to new hypotheses and annotation sharing.

This README provides a package overview, motivation, and installation instructions. For detailed documentation of functions and additional examples, please see the R documentation.

--------------------------------------------------------------------------------
### Installing LOLA

Make sure you have [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html) (bioconductor package) installed:
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```
Then, install the development version directly from github with devtools:
```{r}
devtools::install_github("sheffien/LOLA")
```

Or, clone the repo and install from there:
```{r}
install.packages("path/to/LOLA", repos=NULL)
```
--------------------------------------------------------------------------------
### Running LOLA

For examples and workflows for LOLA, please check out the following [R vignettes](vignettes/) to get you started:

* [Getting Started with LOLA](vignettes/gettingStarted.Rmd)
* [Using the LOLA Core Database](vignettes/usingLOLACore.Rmd) (Requires database files; see below)
* [Choosing a Universe](vignettes/choosingUniverse.Rmd) (Requires database files; see below)

--------------------------------------------------------------------------------
### LOLA Core Database

You can download a core region set database, or (`regionDB`). There are two download options: you can download pre-cached `.RData` files, which LOLA can load in about 30 seconds (requires the [simpleCache R package](http://github.com/sheffien/simpleCache)); or the complete database which additionally includes raw text region files, which LOLA can load and cache in about 30 minutes. LOLA Core currently  contains region sets from hg19 and mm10. We are in the process of adding hg38, among others.

In addition to the LOLA Core database, we also maintain a second database, LOLA Exended, which has additional region sets, which are not as well curated as the Core database (detailed contents are listed below)

The latest LOLA Core and Extended databases can be downloaded here:

* [LOLACore Full database](http://big.databio.org/regionDB/LOLACore_latest.tgz) (Raw source and processed caches, ~1GB)
* [LOLACore Cached database](http://big.databio.org/regionDB/LOLACoreCaches_latest.tgz) (Processed cache files only, ~200Mb)
* [LOLAExtFull database](http://big.databio.org/regionDB/LOLAExt_latest.tgz) (Raw source and processed caches, ~1GB)
* [LOLAExt Cached database](http://big.databio.org/regionDB/LOLAExtCaches_latest.tgz) (Processed cache files only, ~200Mb)
* [Vignette example data](http://big.databio.org/regionDB/lola_vignette_data_150505.tgz) (For testing LOLA Core, ~20Mb)

I recommend using the cached version, unless you need the raw files for something else.
To do this, you'll need to grab [simpleCache](http://github.com/sheffien/simpleCache) (which you may find it useful for other projects, too), also installable with devtools.

```{r}
devtools::install_github("sheffien/simpleCache")
```

Current contents of LOLA core:

* hg19
  1. Transcription Factor binding sites from  [ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/)
  2. Tissue clustered DNase hypersensitive sites from [Sheffield et al. (2013)](http://dnase.genome.duke.edu)
  3. [Codex database](http://codex.stemcells.cam.ac.uk/)
  4. A collection of UCSC feature tables (like CpG islands)
  5. Cistrome database from [Cistrome](http://dx.doi.org/10.1186/gb-2011-12-8-r83)
  6. Epigenome databases from [Cistrome](http://dx.doi.org/10.1186/gb-2011-12-8-r83)
* mm10
  1. [Codex database](http://codex.stemcells.cam.ac.uk/)
  2. Cistrome database (in process)
  3. Encode TFBS (in process)

Current contents of LOLA Extended:
* hg19
  1. Roadmap epigenomics regions
  2. JASPAR motif matches


We're actively adding new collections, so stay tuned. Please contribute! LOLA Core is just the beginning: you can add your own region sets to test enrichment with whatever you like. Here's how to build a custom database:

--------------------------------------------------------------------------------
### Building a custom database

LOLA can read your custom region sets the same way it reads LOLA Core. Check out the raw LOLA Core database for an example of how to organize your own custom database. Your custom database is a bunch of genomic regions, organized first into region sets, and then into collections of region sets.

A bit of terminology:

* Region set: several regions with some shared biological annotation, like a ChIP-seq experiment, represented by a bed file.
* Collection: a named group of bed files

Start by creating collections of bed files: a collection is just a folder with a bunch of bed files in it. For example, you may find a paper that produced 100 data sets you're interested in testing for overlap with your data. To make a collection for the paper, create a folder, name it after the paper, and then put the 100 bed files _into a subfolder called `regions`_. Drop this collection into a parent database folder (perhaps `hg19`) that holds all your collections, and you're good to go!

If you find yourself creating lots of custom collections, you should consider sharing them to improve the LOLA Core database! I'm always looking for additional datasets to add.

##### Basic minimal requirements for a collection

A collection is a folder that contains the following items:

1. `regions/` subfolder with bed-like (chr,start,stop) files inside (**REQUIRED**)
2. `collection.txt` file describing the collection (**RECOMMENDED**)
3. `index.txt` file describing the regions (**RECOMMENDED**)
4. Scripts or descriptions on how to reproduce the collection (**OPTIONAL**)

##### Guidelines for collections

* **All region sets within the collection folder should be in a subfolder named `regions`.**

* For convenience and efficiency, aim for collections between 50 and 2000 region sets. Around 250 is ideal. The software can handle less if you have to, but try lumping small collections together logically, if possible. It will make it easier to organize things in the future. If you lump different sources together, make sure to annotate with appropriate column (see below).

* Name your collection folder something short and informative. The name of the first author of the paper is good, if the collection is completely or mostly derived from a single paper. Otherwise, something general describing all the files.

* Name your bed files something short and informative. If you provide no additional annotation information (see below), this (along with the collection name) will be the only way to identify the region set. No need to put the collection name into the bed name.

##### Annotating collections

You should annotate your collections by putting a file named `collection.txt` into each collection folder. This file should be a 2-line TSV file with a header line including these columns:

* collector (your name)
* date (time you produced the collection)
* source (paper or website where you got the data)
* description (free form field for details)

Example file:

collector		|date		|source		|description
---------------------|-------------|--------------------|-----------
John Doe		|2015-01-03	|Ziller et al.(2014) | Methylation data downloaded from the Ziller paper, files renamed and curated manually.


##### Annotating region sets

You should annotate your region sets by putting a file named `index.txt` into each collection folder. This is not required, but suggested. This file should be a TSV file with a header line including at least a column called `filename`, which points to files in that collection folder. You can then add additional annotation columns; LOLA will recognize and use these column names:

* filename (must match files in the collection exactly)
* description
* cellType
* tissue
* antibody (for ChIP experiments)
* treatment
* dataSource (for publication author, database, etc.)

Any other column names will be ignored, so add whatever else you like. You can also feel free to annotate as many or as few columns, and rows, as you wish. LOLA will simply use as much annotation information as you give it, defaulting to identifying a sample with only the file name if you provide nothing else. So, for example, a `0index` file may look like this:

filename	|cellType	|antibody
--------------|-------------|--------
regionset1.bed|K562		|GATA2
regionset2.bed|K562		|CTCF

These `collection.txt` and `index.txt` annotation files are put inside the collection folder so that a collection is a self-contained entity that can be easily moved.

##### Example custom database

Your folder hierarchy looks something like this:

* regionDB
  * hg19
    * collection1
      * collection.txt
      * index.txt
      * regions/
        * regionset1.bed
        * regionset2.bed
        * regionset3.bed
    * collection2
      * collection.txt
      * index.txt
      * regions/
        * regionset1.bed
        * regionset2.bed
        * regionset3.bed
    * collection3
      * collection.txt
      * regions/
        * bed files...

Then simply pass the `regionDB/hg19` folder (the parent folder containing your collections) to `loadRegionDB()` and it will automatically read and annotate your region collections.

##### Tips
* Your region files really just need the first 3 columns to be chr, start, and end -- no need to follow exact bed specifications.

* Your files don't _have_ to end with `.bed` -- just make sure they are text files. Right now there's no gzip file reading, but this may change in the future.

* You don't have to annotate each file in a collection in the same way, but it's helpful. Just put in whatever you have and LOLA will default to file name for files you don't annotate better.

* You could create your initial `index.txt` file by executing `ls > index.txt` in a collection folder. Now, add a first line containing `filename`, open in spreadsheet software, and start annotating!

* You could stick other annotation files in the parent collection folder if you want. LOLA will ignore them.

* On first load of a collection, LOLA will automatically produce a file called `sizes.txt` containing the size of each set.

* Make sure all files in a collection, and all collections in parent folder, use the same reference genome.

* If you have a single file with different collections (like a segmentation), you can use a function `splitFileIntoCollection()` to divide it into separate bed files so LOLA can understand it.
