## Introduction

This package distributes the code illustrating the use-cases described in the _Using R for proteomics data analysis_ manuscript. 
The package illustrates how R and a selection of dedicated packages can be used to access mass spectrometry base proteomics data, manipulate and visualise it, how to process label-free and labelled quantitative data and how to analyise the quantitation data. 

The package distributes a set of function to access data from the ProteomeXchange `PXD000001` data, used in several examples, as well as a detailed document containing the exact code to reproduce all the analyses presented in the manuscript as well as other application examples. The latter document is called the package vignette, and can be accessed once the package is installed (see below) with the `RforProteomics()` function. Alternatively, the vignette can be downloaded alternatively as a pdf file [here](http://proteome.sysbiol.cam.ac.uk/lgatto/RforProteomics/RforProteomics.pdf). 

## Installation

Since the package is not in Bioconductor, a installation script that will install unmet dependencies and `RforProteomics` automatically is available. From R, type 

```r
source("https://raw.github.com/lgatto/RforProteomics/master/inst/scripts/installR4P.R")
```

## Help

To obtain help or additional information about the `RforProteomics` package, please contact [me](http://proteome.sysbiol.cam.ac.uk/lgatto/). For help about the packages presented in the vignette or manuscript, please refer to the [R mailing list](https://stat.ethz.ch/mailman/listinfo/r-help), [Biodonductor mailing list](http://www.bioconductor.org/help/mailing-list/#bioconductor) (if the packaage is in Bioconductor) and/or the respective package authors. 
