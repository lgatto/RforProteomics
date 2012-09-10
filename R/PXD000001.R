.getPXD000001 <- function(destdir, src, gunzip = TRUE) {
  ans <- FALSE
  on.exit(invisible(ans))
  dest <- basename(src)
  dest <- file.path(destdir, dest)    
  download.file(src, destfile = dest)
  if (gunzip)
    gunzip(dest)
  ans <- TRUE
}


##' Downloads the PXD000001 mzXML file in the \code{destdir}
##' directory. The resulting file is named
##' \code{TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML.gz}.
##'
##' @title Download the PXD000001 mzXML file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns \code{TRUE} if download and gunzip of 
##' the file was successful, \code{FALSE} otherwise. 
##' @author Laurent Gatto
getPXD000001mzXML <- function(destdir = ".") {
  src <- "http://proteome.sysbiol.cam.ac.uk/lgatto/RforProteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML.gz"
  .getPXD000001(destdir, src, TRUE)
}


##' Downloads the PXD000001 mzTba file in the \code{destdir}
##' directory. The resulting file is named
##' \code{F063721.dat-mztab.txt}.
##'
##' @title Download the PXD000001 mzTab file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns \code{TRUE} if download and gunzip of 
##' the file was successful, \code{FALSE} otherwise. 
##' @author Laurent Gatto
getPXD000001mzTab <- function(destdir = ".") {
  src <- "ftp://ftp.pride.ebi.ac.uk/2012/03/PXD000001/F063721.dat-mztab.txt"
  .getPXD000001(destdir, src, FALSE)
}


##' Downloads the PXD000001 mzData file in the \code{destdir}
##' directory. The resulting file is named
##' \code{PRIDE_Exp_Complete_Ac_22134.xml}
##'
##' @title Download the PXD000001 mzTab file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns \code{TRUE} if download and gunzip of 
##' the file was successful, \code{FALSE} otherwise. 
##' @author Laurent Gatto
getPXD000001mzData <- function(destdir = ".") {
  src <- "ftp://ftp.pride.ebi.ac.uk/2012/03/PXD000001/PRIDE_Exp_Complete_Ac_22134.xml.gz"
  .getPXD000001(destdir, src, TRUE)
}



