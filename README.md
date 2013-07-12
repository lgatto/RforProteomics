## Introduction

This package distributes the code illustrating the use-cases described in the 
_Using `R` and Bioconductor for proteomics data analysis_ manuscript
([pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23692960) and [pre-print](http://arxiv.org/abs/1305.6559)). 
The package illustrates how `R` and a selection of dedicated packages can be used 
to access mass-spectrometry proteomics data, manipulate and visualise it, 
how to process label-free and labelled quantitative data and how to analyse the quantitation data. 

The package will be updated beyond the content of the manuscript to keep up-to-date with progress in the area.
The [github page](https://github.com/lgatto/RforProteomics) can be used to edit the wiki and file new issues 
related to the package itself of general needed for proteomics that should be addressed in `R`.

It would be great if this work could stimulate a wider participation to use `R` and develop `R`
packages for proteomics and promote interaction between computational biologists working in the
field of proteomics, in particular by facilitating interoperability between their software.
The [rbioc-sig-proteomics](https://groups.google.com/forum/#!forum/rbioc-sig-proteomics) group
has tentatively been set up to provide an question and discussion forum for interested parties.
You can also follow updates to some of the relevant packages on twitter via [@RforProteomics](https://twitter.com/r4proteomics).
Do not hesitate to get in touch with further suggestions.


## Data and vignette

The package also distributes a set of function to access data from the [ProteomeXchange](http://www.proteomexchange.org/) [`PXD000001`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000001) data, used in several examples, as well as a detailed document containing the exact code to reproduce all the analyses presented in the manuscript as well as other application examples. The latter document is called the package vignette, and can be accessed once the package is installed (see below) with the `RforProteomics()` function. Alternatively, the vignette can be downloaded as a pdf file [here](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.pdf). 

## Installation

The package is now available in [Bioconductor](http://bioconductor.org/packages/devel/data/experiment/html/RforProteomics.html) version 2.13. To install the package and its documentation, start `R` (>= `3.0.0`) and enter:

```c
source("http://bioconductor.org/biocLite.R")
biocLite("RforProteomics")
```

To install all dependencies (78 packages, including `RforProteomics`) and reproduce the code in the vignette, replace the last line in the code chunk above with:

```c
biocLite("RforProteomics", dependencies = TRUE)
```

Some packages might not be available in the stable Bioconductor version (currently 2.12).
To have access to new packages that are available in the development branch or make use of very recent features of packages in the stable version, one would need to use the development version of Bioconductor by setting the `useDevel` variable as shown below (or use the `R` development version, only recommended for advanced users though).

```c
library("BiocInstaller")
useDevel(TRUE)
```

## Help

To obtain help or additional information about the `RforProteomics` package, please contact [me](http://proteome.sysbiol.cam.ac.uk/lgatto/). For help about the packages presented in the vignette or manuscript, please refer to the [R mailing list](https://stat.ethz.ch/mailman/listinfo/r-help), [Bioconductor mailing list](http://www.bioconductor.org/help/mailing-list/#bioconductor) (if the package is in Bioconductor) and/or the respective package authors. 

For general resources about `R`, see the corresponding section in the vignettes and the [TeachingMaterial](https://github.com/lgatto/TeachingMaterial) repository.
