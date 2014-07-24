##' Opens the package vignettes.
##'
##' @title Opens RforProteomics vignettes
##' @return An instance of class \code{vignette},
##' as returned by \code{\link{vignette}}
##' @author Laurent Gatto
RforProteomics <- function()
  vignette("RforProteomics", package = "RforProteomics")



##' Finds the package names that have a specific biocView.
##'
##' @title Packages in a biocView
##' @param view The biocView of interest. For example
##' \code{"Proteomics"}.
##' @param rep Repository of interest. One of \code{"BioCsoft"},
##' \code{"BioCann"}, \code{"BioCexp"} or \code{"BioCextra"}.
##' @param biocVersion A \code{character} with the Bioconductor
##' version of interest. For example \code{"2.14"}.
##' @return An instance of class \code{BiocView}. \code{NULL} if the
##' the biocView was not found.
##' @author Laurent Gatto
getPackagesInBiocView <- function(view, 
                                  rep = c("BioCsoft", "BioCann",
                                      "BioCexp", "BioCextra"),
                                  biocVersion) {
    suppressPackageStartupMessages(require("biocViews"))
    data(biocViewsVocab)
    rep <- match.arg(rep)
    biocMirror <- getOption("BioC_mirror", "http://bioconductor.org")
    biocPaths <- switch(rep,         
                        BioCsoft = "bioc",
                        BioCann = "data/annotation", 
                        BioCexp = "data/experiment",
                        BioCextra = "extra")   
    rep <- paste(biocMirror,
                 "packages",
                 biocVersion,
                 biocPaths, 
                 sep = "/") 
    bv <- getBiocViews(rep, biocViewsVocab, "NoViewProvided")
    if (!view %in% names(bv)) {
        message("BiocView ", view, " not found.")
        return(NULL)
    }  
    return(bv[[view]])
}


##' Format a \code{BiocView} as a \code{data.frame}.
##'
##' @title Package descriptions
##' @param x An instance of class \code{BiocView}, as produced by
##' \code{getPackagesInBiocView}.
##' @param nsub A \code{logical} indicating \code{"\n"} are to be
##' replaced by a space.
##' @param version A \code{logical} specifying if the package version
##' should be added.
##' @return A \code{data.frame} with package information. 
##' @author Laurent Gatto
packageDF <- function(x, nsub = TRUE, version = TRUE) {
  Package <- sapply(x@packageList, function(x) x@Package)
  Title <- sapply(x@packageList, function(x) x@Title)
  Version <- sapply(x@packageList, function(x) x@Version)
  if (nsub)
    Title <- sub("\n", " ", Title)
  ans <- data.frame(Package, Title)
  if (version) ans$Version <- Version
  return(ans)
}



proteomicsPackages <- function(biocv) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    packageDF(getPackagesInBiocView("Proteomics",
                                    biocVersion = biocv))
}

massSpectrometryPackages <- function(biocv) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    packageDF(getPackagesInBiocView("MassSpectrometry",
                                    biocVersion = biocv))
}

massSpectrometryDataPackages <- function(biocv) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    packageDF(getPackagesInBiocView("MassSpectrometryData",
                                    rep = "BioCexp",
                                    biocVersion = biocv))
}
