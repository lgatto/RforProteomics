\name{getPXD000001mzXML}
\alias{getPXD000001mzXML}
\title{Download the PXD000001 mzXML file}
\usage{
  getPXD000001mzXML(destdir = ".")
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
  Downloads the PXD000001 mzXML file in the \code{destdir}
  directory. The resulting file is named
  \code{TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML.gz}.
}
\author{
  Laurent Gatto
}

