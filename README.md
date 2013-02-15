## Introduction

This package distributes the code illustrating the use-cases described in the 
_Using `R` and Bioconductor for proteomics data analysis_ manuscript. 
The package illustrates how `R` and a selection of dedicated packages can be used 
to access mass spectrometry proteomics data, manipulate and visualise it, 
how to process label-free and labelled quantitative data and how to analyse the quantitation data. 

The package will be updated beyond the content of the manuscript to keep up-to-date with progress in the area.

## Data and vignette

The package also distributes a set of function to access data from the [ProteomeXchange](http://www.proteomexchange.org/) [`PXD000001`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000001) data, used in several examples, as well as a detailed document containing the exact code to reproduce all the analyses presented in the manuscript as well as other application examples. The latter document is called the package vignette, and can be accessed once the package is installed (see below) with the `RforProteomics()` function. Alternatively, the vignette can be downloaded as a pdf file [here](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.pdf). 

## Installation

The package is now available in [Bioconductor](http://bioconductor.org/packages/devel/data/experiment/html/RforProteomics.html) version 2.12. To install the package, start `R` (currently devel version only) and enter:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("RforProteomics")
```

As of writing, Bioc 2.12 is the development branch, which require the development version of `R` to install packages using `biocLite`. If you have an earlier version of `R` (`R-2.15.2`, the latest stable version is recommended), start by installing the package dependencies as shown below. The `deps` variable is a list of all packages that are needed to replicate all the code illustrated in the package. The second line loads the `installPackages` function, available on a server. This function will first check if any of the packages listed in `deps` are not already available and proceed only with missing dependencies. The installation uses the recommended `BiocInstaller` package, which is installed if not yet available. If packages are out dated, they will be updated.

```r
deps <- c('R.utils', 'Biobase', 'mzR', 'MSnbase', 
          'xcms', 'msdata', 'isobar', 'MALDIquant', 'MALDIquantForeign',
          'readBrukerFlexData', 'synapter', 'synapterdata', 
          'IPPD', 'Rdisop', 'OrgMassSpecR', 'BRAIN', 
          'rols', 'hpar', 'GO.db', 'org.Hs.eg.db', 
          'biomaRt', 'RColorBrewer', 'ggplot2', 'reshape2', 
          'knitr')
source("http://proteome.sysbiol.cam.ac.uk/lgatto/src/installPackages.R")
installPackages(deps)
```

Download the appropriate package from the [Bioconductor landing page](http://bioconductor.org/packages/devel/data/experiment/html/RforProteomics.html) and install manually using `install.packages(..., repos = "NULL")` or the GUI front-end of your favourite R interface.

Once installed, the package with loaded with

```r
library("RforProteomics")
```

If you are using `R-2.15.x`, a message will warn that `RforProteomics` has been built using `R-3.0.0` (the development version). 

## Help

To obtain help or additional information about the `RforProteomics` package, please contact [me](http://proteome.sysbiol.cam.ac.uk/lgatto/). For help about the packages presented in the vignette or manuscript, please refer to the [R mailing list](https://stat.ethz.ch/mailman/listinfo/r-help), [Bioconductor mailing list](http://www.bioconductor.org/help/mailing-list/#bioconductor) (if the package is in Bioconductor) and/or the respective package authors. 
