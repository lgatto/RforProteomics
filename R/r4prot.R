##' Opens the package vignettes.
##'
##' @title Opens the RforProteomics vignette
##' @return An instance of class \code{vignette}. Used for its side
##' effect, opening the vignette.
##' @author Laurent Gatto
RforProteomics <- function()
  vignette("RforProteomics", package = "RforProteomics")

##' Opens the visualisation vignette
##'
##' @title Opens the visualisation vignette
##' @return An instance of class \code{vignette}. Used for its side
##' effect, opening the vignette.
##' @aliases RProtViz
##' @author Laurent Gatto
RProtVis <- function()
    vignette("RProtVis", package = "RforProteomics")

RProtViz <- RProtVis

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
    biocViewsVocab <- NULL ## no visible binding warning
    data(biocViewsVocab, package = "biocViews", envir = environment())
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



proteomicsPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    if (cache) {
        f <- dir(system.file("extdata", package = "RforProteomics"),
                 full.names = TRUE, pattern = "lpp.rds")
        x <- readRDS(f)
        x <- x[[biocv]]
        if (is.null(x))
            warning(biocv, " not cached, retrieving from Bioconductor.org",
                    call. = FALSE, immediate. = TRUE)
        else return(x)
    }
    packageDF(getPackagesInBiocView("Proteomics",
                                    biocVersion = biocv))
}

massSpectrometryPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    if (cache) {
        f <- dir(system.file("extdata", package = "RforProteomics"),
                 full.names = TRUE, pattern = "lmsp.rds")
        x <- readRDS(f)
        x <- x[[biocv]]
        if (is.null(x))
            warning(biocv, " not cached, retrieving from Bioconductor.org",
                    call. = FALSE, immediate. = TRUE)
        else return(x)
    }
    packageDF(getPackagesInBiocView("MassSpectrometry",
                                    biocVersion = biocv))
}

massSpectrometryDataPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocInstaller::biocVersion())
    if (cache) {
        f <- dir(system.file("extdata", package = "RforProteomics"),
                 full.names = TRUE, pattern = "lmsdp.rds")
        x <- readRDS(f)
        x <- x[[biocv]]
        if (is.null(x))
            warning(biocv, " not cached, retrieving from Bioconductor.org",
                    call. = FALSE, immediate. = TRUE)
        else return(x)
    }
    packageDF(getPackagesInBiocView("MassSpectrometryData",
                                    rep = "BioCexp",
                                    biocVersion = biocv))
}


msDataTab <- function() {
    dat <- c("Raw", "Identification", "Identification", "Quantitative", "Raw", "Peak lists")
    frmt <- c("mzXML or mzML", "mzIdentML", "mzIdentML", "mzTab", "mzML or mzXML", "mgf")
    robj <- c("mzRpwiz and mzRramp", "mzRident", "mzID", "MSnSet", "MSnExp", "MSnExp")
    pkg <- c("mzR", "mzR", "mzID", "MSnbase", "MSnbase", "MSnbase")
    ans <- data.frame(dat, frmt, robj, pkg)
    colnames(ans) <- c("Data type", "File format", "Data structure", "Package")
    ans
}
