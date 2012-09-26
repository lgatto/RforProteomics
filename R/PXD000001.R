.getPXD000001 <- function(destdir, src, unpack = TRUE) {
  dest <- basename(src)
  dest <- file.path(destdir, dest)
  dest2 <- gsub("[.]gz$", "", dest)
  if (!file.exists(dest2)) {
    if (!file.exists(dest))
      download.file(src, destfile = dest)
    if (unpack)
      gunzip(dest)
  }
  invisible(dest2)
}


##' Unless already present, downloads the PXD000001 mzXML file
##' in the \code{destdir} directory. The resulting file is named
##' \code{TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML}.
##'
##' @title Download the PXD000001 mzXML file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns the name of the downloaded file.
##' @author Laurent Gatto
getPXD000001mzXML <- function(destdir = ".") {
  src <- "http://proteome.sysbiol.cam.ac.uk/lgatto/RforProteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML.gz"
  .getPXD000001(destdir, src, TRUE)
}


##' Unless already present, downloads the PXD000001 mzTab file
##' n the \code{destdir} directory. The resulting file is named
##' \code{F063721.dat-mztab.txt}.
##'
##' @title Download the PXD000001 mzTab file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns the name of the downloaded file.
##' @author Laurent Gatto
getPXD000001mzTab <- function(destdir = ".") {
  src <- "ftp://ftp.pride.ebi.ac.uk/2012/03/PXD000001/F063721.dat-mztab.txt"
  .getPXD000001(destdir, src, FALSE)
}


##' Unless already present, downloads the PXD000001 mzData file
##' in the \code{destdir} directory. The resulting file is named
##' \code{PRIDE_Exp_Complete_Ac_22134.xml}
##'
##' @title Download the PXD000001 mzTab file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns the name of the downloaded file.
##' @author Laurent Gatto
getPXD000001mzData <- function(destdir = ".") {
  src <- "ftp://ftp.pride.ebi.ac.uk/2012/03/PXD000001/PRIDE_Exp_Complete_Ac_22134.xml.gz"
  .getPXD000001(destdir, src, TRUE)
}



