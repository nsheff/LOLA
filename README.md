LOLA: Location Overlap Analysis
-------------------------------

LOLA is an R package providing functions for testing overlap of sets of genomic regions with public and custom databases. This README provides a package overview, motivation, and installation instructions. For detailed documentation of functions and additional examples, please see the R documentation.

LOLA is deployable anywhere, using generic custom bed-file database.

--------------------------------------------------------------------------------
### Installing LOLA
Install the development version directly from github with devtools:
```{r}
require(devtools)
install_github("sheffien/LOLA") 
```

Or, clone the repo and install from there:
```{r}
packageFolder = "~/R/LOLA";
install.packages(packageFolder, repos=NULL)
```

--------------------------------------------------------------------------------
### Running LOLA

For examples and workflows for `LOLA`, please check out the [R vignettes](vignettes/) to get you started:

* [Getting Started with LOLA](vignettes/gettingStarted.Rmd)

--------------------------------------------------------------------------------
### LOLA Core Database

LOLA comes with a core database (`regionDB`) that includes ENCODE transcription factor binding sites, the cistrome database, and DNase hypersensitive sites. There are two download options: you can download pre-cached `.RData` files, which `LOLA` can load in about 30 seconds (requires the [simpleCache R package](http://github.com/sheffien/simpleCache)); or the complete database which additionally includes raw text region files, which `LOLA` can load in about 30 minutes. LOLA Core currently **only contains region sets from hg19, but we will be adding mm10 at some point**.

The latest LOLA Core database can be downloaded here:

* [Full database](http://www.biomedical-sequencing.at/bocklab/nsheffield/regionDB/regionDB_150416.tgz) (Raw source and processed caches - 840Mb)
* [Cached database](http://www.biomedical-sequencing.at/bocklab/nsheffield/regionDB/regionDBcaches_150416.tgz) (Processed cache files only - 167Mb)

If you're just using `LOLA`, I recommend using the cached version, unless you need the raw files for something else. 
To do this, you'll need to grab my R package `simpleCache` (which you may find it useful for other projects, too).

```{r}
install_github("sheffien/simpleCache")
```

`LOLA` Core is just the beginning: you can add your own region sets to test enrichment with whatever you like. Here's how to build a custom database:

--------------------------------------------------------------------------------
### Building a custom database

LOLA provides functions that can read your bed files, building a custom database. Start by creating collections of bed files: a collection is just a folder with a bunch of bed files in it. For example, you may find a paper that produced 100 data sets you're interested in testing for overlap with your data. To make a collection for the paper, create a folder, name it after the paper, and then put the 100 bed files _into a subfolder called `regions`_. Drop this collection into a parent database folder that holds all your collections, and you're good to go!

##### Guidelines for collections

* **All regions within the collection folder should be in a subfolder named `regions`.**

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
* cell-type
* tissue
* antibody (for ChIP experiments)
* treatment
* data-source (for publication author, database, etc.)

Any other column names will be ignored, so add whatever else you like. You can also feel free to annotate as many or as few columns, and rows, as you wish. LOLA will simply use as much annotation information as you give it, defaulting to identifying a sample with only the file name if you provide nothing else. So, for example, a `0index` file may look like this:

filename	|cell-type	|antibody
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
