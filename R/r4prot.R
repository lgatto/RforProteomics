##' Opens the package and visualisation vignettes.
##'
##' @title Opens the package vignettes
##' @return An instance of class `vignette`. Used for its side effect,
##'     opening a vignette.
##' @author Laurent Gatto
##' @md
##' @export
RforProteomics <- function()
  vignette("RforProteomics", package = "RforProteomics")

##' @export
##' @rdname RforProteomics
RProtVis <- function()
    vignette("RProtVis", package = "RforProteomics")

##' @export
##' @rdname RforProteomics
RProtViz <- RProtVis

##' Finds the package names that have a specific biocView.
##'
##' @title Packages in a biocView
##' @param view The biocView of interest. For example `"Proteomics"`.
##' @param rep Repository of interest. One of `"BioCsoft"`,
##'     `"BioCann"`, `"BioCexp"` or `"BioCextra"1.
##' @param biocVersion A `character` with the Bioconductor version of
##'     interest. For example `"2.14"`.
##' @return An instance of class `BiocView`. `NULL` if the the
##'     biocView was not found.
##' @md
##' @author Laurent Gatto
##' @import biocViews getBiocViews
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



##'   Searches for all the packages with the \code{"Proteomics"}
##'   (software), \code{"MassSpectrometry"} (software) and
##'   \code{"MassSepctrometryData"} (data) packages and return their
##'   names, titles and versions as a \code{data.frame}. The
##'   (unexported but documented) underlying functions are
##'   \code{RforProteomics:::getPackagesInBiocView} (to find relevant
##'   package) and \code{RforProteomics:::packageDF}
##'   (\code{data.frame} formatting).
##'
##' @title Proteomics and MS biocView packages
##' @param biocv A \code{character} with the Bioconductor version to
##'     search for relevant packages. If missing, the running version
##'     is used.
##' @param cache A \code{logical} indicating whether cached package
##'     information should be used. Default is \code{FALSE}. All
##'     except development versions are up-to-date.
##' @return A \code{data.frame} with the respective package names,
##'     titles and versions.
##' @author Laurent Gatto
##' @export
##' @examples
##' head(pp <- proteomicsPackages("3.0"))
##' ppc <- proteomicsPackages("3.0", cache = TRUE)
##' all.equal(pp, ppc)
proteomicsPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocManager::version())
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

##' @rdname proteomicsPackages
##' @export
massSpectrometryPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocManager::version())
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

##' @rdname proteomicsPackages
##' @export
massSpectrometryDataPackages <- function(biocv, cache=FALSE) {
    if (missing(biocv))
        biocv <- as.character(BiocManager::version())
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
    dat <- c("Raw", "Raw", "Raw", "Identification", "Identification",
             "Quantitative", "Peak lists", "Imaging", "Imaging")
    frmt <- c("mzXML or mzML", "mzXML or mzML", "mzXML or mzML",
              "mzIdentML", "mzIdentML", "mzTab", "mgf",
              "imzML or Analyze 7.5", "imzML or Analyze 7.5")
    robj <- c("mzRpwiz or mzRramp", "list of MassSpectrum objects",
              "MSnExp", "mzRident", "mzID", "MSnSet", "MSnExp",
              "MSImageSet",  "list of MassSpectrum objects")
    pkg <- c("mzR", "MALDIquantForeign", "MSnbase", "mzR", "mzID",
             "MSnbase", "MSnbase", "Cardinal", "MALDIquantForeign")
    ans <- data.frame(dat, frmt, robj, pkg)
    colnames(ans) <- c("Data type", "File format", "Data structure", "Package")
    ans
}
