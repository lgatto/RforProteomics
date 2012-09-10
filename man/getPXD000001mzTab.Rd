\name{getPXD000001mzTab}
\alias{getPXD000001mzTab}
\title{Download the PXD000001 mzTab file}
\usage{
  getPXD000001mzTab(destdir = ".")
}
\arguments{
  \item{destdir}{A \code{character} with the destination
  folder.}
}
\value{
  Invisibly returns \code{TRUE} if download and gunzip of
  the file was successful, \code{FALSE} otherwise.
}
\description{
  Downloads the PXD000001 mzTba file in the \code{destdir}
  directory. The resulting file is named
  \code{F063721.dat-mztab.txt}.
}
\author{
  Laurent Gatto
}

