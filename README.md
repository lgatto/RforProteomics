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

It would be great if this work could stimulate a wider participation
to use `R` and develop `R` packages for proteomics and promote
interaction between computational biologists working in the field of
proteomics, in particular by facilitating interoperability between
their software.  The
[rbioc-sig-proteomics](https://groups.google.com/forum/#!forum/rbioc-sig-proteomics)
group has tentatively been set up to provide a forum for questions and
discussion for interested parties. Do not hesitate to
[get in touch](http://proteome.sysbiol.cam.ac.uk/lgatto/) for
questions, comments or further suggestions. Note taking about
plans/ideas/direction for R/Bioc and proteomics can be contributed to
the
[`RforProteomics` wiki](https://github.com/lgatto/RforProteomics/wiki).

## Data and vignette

The package uses the dataset
[`PXD000001`](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000001)
from the [ProteomeXchange](http://www.proteomexchange.org/) repository
in several examples. The data can be queries and downloaded from `R`
with the
[`rpx`](http://bioconductor.org/packages/release/bioc/html/rpx.html)
package. The `RforProteomics` vignette is a detailed document
containing the exact code to reproduce all the analyses presented in
the manuscript as well as other application examples. It can be
accessed once the package is installed (see below) with the
`RforProteomics()` function. Alternatively, the vignette can be
downloaded as a pdf file
[here](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.pdf).

A second vignette, `RProtVis` focuses on the visualisation of mass
spectrometry and proteomics data with `R` and Bioconductor. From `R`,
it is currently only available with Bioconductor `>= 3.0` using the
`RProtViz()` function. It can also be consulted
[on-line](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/RProtVis.html)
on the 'RforProteomics' development version page.

## Installation

The package is available on
[Bioconductor](http://bioconductor.org/packages/devel/data/experiment/html/RforProteomics.html)
(version >= 2.13). To install the package and its documentation, start
`R` (>= `3.0.0` required) and type:

```c
source("http://bioconductor.org/biocLite.R")
biocLite("RforProteomics")
```

To install all dependencies (75+ packages, including `RforProteomics`)
and fully reproduce the code in the vignettes, replace the last line
in the code chunk above with:

```c
biocLite("RforProteomics", dependencies = TRUE)
```

## Help

To obtain help or additional information about the `RforProteomics`
package, please contact
[me](http://proteome.sysbiol.cam.ac.uk/lgatto/). For help about the
packages presented in the vignette or manuscript, please refer to the
[R mailing list](https://stat.ethz.ch/mailman/listinfo/r-help),
[Bioconductor mailing list](http://www.bioconductor.org/help/mailing-list/#bioconductor)
(if the package is in Bioconductor) and/or the respective package
authors.

For general resources about `R`, see the corresponding section in the
vignettes and the
[TeachingMaterial](https://github.com/lgatto/TeachingMaterial)
repository.
