\name{proteomicsPackages}
\alias{proteomicsPackages}
\alias{massSpectrometryPackages}
\alias{massSpectrometryDataPackages}


\title{
  Proteomics and MS biocView packages
}

\description{

  Searches for all the packages with the \code{"Proteomics"} (software),
  \code{"MassSpectrometry"} (software) and \code{"MassSepctrometryData"}
  (data) packages and return their names, titles and versions as a
  \code{data.frame}. The (unexported but documented) underlying
  functions are \code{RforProteomics:::getPackagesInBiocView} (to find
  relevant package) and \code{RforProteomics:::packageDF}
  (\code{data.frame} formatting).  }

\usage{
proteomicsPackages(biocv, cache=FALSE)
massSpectrometryPackages(biocv, cache=FALSE)
massSpectrometryDataPackages(biocv, cache=FALSE)
}

\arguments{
  \item{biocv}{A \code{character} with the Bioconductor version to
    search for relevant packages. If missing, the running version is used.
  }
  \item{cache}{A \code{logical} indicating whether cached package
  information should be used. Default is \code{FALSE}. All except
  development versions are up-to-date.}
}

\value{
  A \code{data.frame} with the respective package names, titles and
  versions. 
}

\author{
  Laurent Gatto
}

\examples{
head(pp <- proteomicsPackages("3.0"))
ppc <- proteomicsPackages("3.0", cache = TRUE)
all.equal(pp, ppc)
}
