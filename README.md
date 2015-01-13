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
### Building a custom database

LOLA provides functions that can read your bed files. Just create a folder for your database, then put collections of bed files into subfolders. Actually, your files just need the first 3 columns to be chr, start, and end -- no need to follow exact bed specifications. So your folder hierarchy looks something like this:

* regionDB
  * hg19
    * collection1
      * 0index
      * file1.bed
      * file2.bed
      * file3.bed
    * collection2
    * collection 3

Then simply pass the hg19 folder to loadRegionDB() and it will automatically read and annotate your region collections.

You can (optionally) annotate your regions by putting a file named `0index` into a collection folder. This file should be a TSV file with a header line including at least a column called `filename`, which points to files in that collection folder. You can then add additional annotation columns; LOLA will recognize and use these column names:

* filename
* description
* source 
* antibody
* treatment

Any other column names will be ignored, so add whatever else you like. You can also feel free to annotate as many or as few columns, and rows, as you wish. LOLA will simply use as much annotation information as you give it, defaulting to identifying a sample with only the file name if you provide nothing else. So, for the above structure, a 0index may look like this:

filename	|source	|antibody
file1.bed	|K562		|GATA2 
file3.bed	|K562		|CTCF



