## Introduction

The
[`RforProteomics` package](http://www.bioconductor.org/packages/release/data/experiment/html/RforProteomics.html)
distributes and extends the use-cases described in the _Using `R` and
Bioconductor for proteomics data analysis_ manuscript
([pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23692960) and
[pre-print](http://arxiv.org/abs/1305.6559)).  The package illustrates
how `R` and a selection of dedicated packages that can be used to
access mass-spectrometry proteomics data, manipulate and visualise it,
how to process label-free and labelled quantitative data and how to
analyse the quantitation data.

The package will be updated beyond the content of the manuscript to
keep up-to-date with progress in the area.  The
[github page](https://github.com/lgatto/RforProteomics) can be used to
edit the wiki and file new issues related to the package itself of
general needed for proteomics that should be addressed in `R`.

> **NB**: I you are interested in R packages for mass
> spectrometry-based proteomics and metabolomics, see also the [R for
> Mass Spectrometry initiative](https://www.rformassspectrometry.org/)
> packages and the [tutorial
> book](https://rformassspectrometry.github.io/docs/)


It would be great if this work could stimulate a wider participation
to use `R` and develop `R` packages for proteomics and promote
interaction between computational biologists working in the field of
proteomics, in particular by facilitating interoperability between
their software.  The official [Bioconductor support
site](https://support.bioconductor.org/) is the channel of choice to
ask questions about specific Bioconductor packages.

## Data and vignette

The package uses the
dataset
[`PXD000001`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000001) from
the [ProteomeXchange](http://www.proteomexchange.org/) repository in
several examples. The data can be queried and downloaded from `R` with
the
[`rpx`](http://bioconductor.org/packages/release/bioc/html/rpx.html)
package. The `RforProteomics` vignette is a detailed document
containing the exact code to reproduce all the analyses presented in
the manuscript as well as other application examples. It can be
accessed once the package is installed (see below) with the
`RforProteomics()` function. Alternatively, the vignettes can be read
online
[here](https://lgatto.github.io/RforProteomics/articles/RforProteomics.html) and
[here](https://lgatto.github.io/RforProteomics/articles/RProtVis.html).

A second vignette, `RProtVis` focuses on the visualisation of mass
spectrometry and proteomics data with `R` and Bioconductor. From `R`,
it is currently only available with Bioconductor `>= 3.0` using the
`RProtViz()` function. It can also be consulted
[on-line](http://bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RProtVis.html)
on the 'RforProteomics' development version page.

## Installation

The package is available on
[Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/RforProteomics.html)
(version >= 2.13). To install the package and its documentation, start
`R` (>= `3.0.0` required) and type:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("RforProteomics")
```

To install all dependencies (75+ packages, including `RforProteomics`)
and fully reproduce the code in the vignettes, replace the last line
in the code chunk above with:

```
BiocManager::install("RforProteomics", dependencies = TRUE)
```

## Collaborative editing

The community and package authors are invited to contribute to the
package. If you have or know of a package of interest, please
[fork](https://help.github.com/articles/fork-a-repo) the
[repository](https://github.com/lgatto/RforProteomics), add a new
section to the vignette and send a
[pull request](https://help.github.com/articles/creating-a-pull-request). If
you update the vignette, please also add yourself as a contributor to
the package.

## Help

To obtain help or additional information about the `RforProteomics`
package and about the packages presented in the vignette or
manuscript, please use the [Bioconductor support
site](https://support.bioconductor.org/).

For general resources about `R`, see the corresponding section in the
vignettes and the
[TeachingMaterial](https://github.com/lgatto/TeachingMaterial)
repository.
