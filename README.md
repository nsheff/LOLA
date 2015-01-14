LOLA: Location Overlap Analysis
-------------------------------

LOLA is an R package providing functions for testing overlap of a set of genomic regions with public and custom databases. This REAMDE provides a package overview, installation instructions, and use cases. For detailed documentation of functions and additional examples, please see the R documentation.

--------------------------------------------------------------------------------
### Installing LOLA
install the development version directly from github with devtools

To install LOLA with devtools (LOLA relies on simpleCache):

```
require(devtools)
install_github("sheffien/simpleCache")
install_github("sheffien/LOLA") 
```

To install a local copy:

```
packageFolder = "~/Rncs/LOLA";
install.packages(packageFolder, repos=NULL)
```

--------------------------------------------------------------------------------
### Running LOLA

```
library(LOLA)
?LOLA
```

Now take a look at examples.R for some examples. Vignettes are forthcoming.

--------------------------------------------------------------------------------
### LOLA Core

LOLA comes with a core database that includes encode transcription factor binding sites, the cistrome database, and DNase hypersensitive sites. You may want to add additional region sets specific to your application (instructions follow).

--------------------------------------------------------------------------------
### Building a custom database

LOLA provides functions that can read your bed files, building a custom database. Start by creating collections of bed files: a collection is just a folder with a bunch of bed files in it. For example, you may find a paper that produced 100 data sets you're interested in testing for overlap with your data. To make a collection for the paper, create a folder, name it after the paper, and then put the 100 bed files into the folder. Drop this collection into a parent database folder that holds all your collections, and you're good to go!

##### Guidelines for collections

* For convenience and efficiency, aim for collections between 50 and 2000 region sets. Around 250 is ideal. The software can handle less if you have to, but try lumping small collections together logically, if possible. It will make it easier to organize things in the future. If you lump different sources together, make sure to annotate with appropriate column (see below).

* Name your collection something short and informative. The name of the first author of the paper is good, if the collection is completely or mostly derived from a single paper. Otherwise, something general describing all the files.

* Name your bed files something short and informative. If you provide no additional annotation information (see below), this (along with the collection name) will be the only way to identify the region set. No need to put the collection name into the bed name.

##### Annotating regions

You should annotate your regions by putting a file named `0index` into each collection folder. This is not required, but suggested. This file should be a TSV file with a header line including at least a column called `filename`, which points to files in that collection folder. You can then add additional annotation columns; LOLA will recognize and use these column names:

* filename (must match files in the collection exactly)
* description
* cell-type
* tissue
* antibody (for ChIP experiments)
* treatment
* data-source (for publication author, database, etc.)

Any other column names will be ignored, so add whatever else you like. You can also feel free to annotate as many or as few columns, and rows, as you wish. LOLA will simply use as much annotation information as you give it, defaulting to identifying a sample with only the file name if you provide nothing else. So, for the above structure, a `0index` file may look like this:

filename	|source	|antibody
--------------|-------------|--------
file1.bed	|K562		|GATA2 
file3.bed	|K562		|CTCF

The `0` in `0index` simply makes this file sort at the top of the list so you can find it easily. It's located inside the collection folder so that a collection is a self-contained entity that can be easily moved.

##### Example custom database

Your folder hierarchy looks something like this:

* regionDB
  * hg19
    * collection1
      * 0index
      * file1.bed
      * file2.bed
      * file3.bed
    * collection2
      * bed files...
    * collection3
      * bed files...

Then simply pass the `regionDB/hg19` folder (the parent folder containing your collections) to `loadRegionDB()` and it will automatically read and annotate your region collections.

##### Tips
* Your files really just need the first 3 columns to be chr, start, and end -- no need to follow exact bed specifications.

* You don't have to annotate each file in a collection in the same way, but it's helpful. Just put in whatever you have and LOLA will default to file name for files you don't annotate better.

* You could create your initial `0index` file by executing `ls > 0index` in a collection folder. Now, add a first line containing `filename`, open in spreadsheet software, and start annotating!

* Any file starting with a `0` will be ignored by LOLA. So don't name your bed files like this!

* On first load of a collection, LOLA will automatically produce a file called `0sizes` containing the size of each set.
