\name{getPXD000001mzData}
\alias{getPXD000001mzData}
\title{Download the PXD000001 mzTab file}
\usage{
  getPXD000001mzData(destdir = ".")
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
  Downloads the PXD000001 mzData file in the \code{destdir}
  directory. The resulting file is named
  \code{PRIDE_Exp_Complete_Ac_22134.xml}
}
\author{
  Laurent Gatto
}

