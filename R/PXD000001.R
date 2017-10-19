##' Unless already present, downloads \code{src} in the \code{destdir}
##' directory.
##'
##' @title Download a file
##' @param src The url of the file to download.
##' @param destdir The destination directory. Default is \code{"."}.
##' @param unpack Should \code{src} be uncompressed? Default is
##' \code{TRUE}.
##' @param ... Additional paramters passed to
##' \code{\link{download.file}}.
##' @return Invisible returns the full path of the downloaded file.
##' @author Laurent Gatto
downloadData <- function(src, destdir = ".", unpack = TRUE, ...) {
  dest <- basename(src)
  dest <- file.path(destdir, dest)
  dest2 <- gsub("[.]gz$", "", dest)
  if (!file.exists(dest2)) {
    if (!file.exists(dest))
      download.file(src, destfile = dest, ...)
    if (unpack)
      gunzip(dest)
  }
  invisible(dest2)
}

##' Downloads on of multiple Thermo Hela/PRTC data files.
##' 
##' @title Dowload Thermo Hela PRTC data
##' @param src The name of the file to be downloaded. If missing, a
##' vector of possible filenames is returned. If \code{"all"}, all
##' files are downloaded. Alternatively, a pattern can be used to
##' \code{grep} the files from the output \code{getThermoHelaPRTC()} the
##' files to be downloaded.
##' @param destdir Destination directory. Default is \code{"."}.
##' @return Invisibly return the path of the downloaded files. 
##' @author Laurent Gatto
##' @seealso \code{downloadData}
##' @examples
##' getThermoHelaPRTC()
##' getThermoHelaPRTC("design")
##' \dontrun{
##' getThermoHelaPRTC("all")
##' }
getThermoHelaPRTC <- function(src, destdir = ".") {
    url <- "http://proteome.sysbiol.cam.ac.uk/lgatto/RforProteomics/"
    Thermo_Hela_files <-c(
        "design.txt",
        "swissprot_human_canonical_19_09_12.fasta.gz", 
        "Thermo_Hela_PRTC_1.mgf.gz", 
        "Thermo_Hela_PRTC_2.mgf.gz", 
        "Thermo_Hela_PRTC_3.mgf.gz")
    if (missing(src)) {
        return(Thermo_Hela_files)
    } else {
        Thermo_Hela_files2 <-
            paste(url, Thermo_Hela_files, sep = "/")
        if (src == "all") {
            ans <- sapply(Thermo_Hela_files2[1], downloadData, unpack = FALSE)
            ans <- sapply(Thermo_Hela_files2[2:5], downloadData, unpack = TRUE)
        } else {
            Thermo_Hela_files2 <- 
                paste(url, match.arg(src, Thermo_Hela_files), sep = "/")
            ans <- mapply(downloadData,
                          Thermo_Hela_files2,
                          unpack = grepl("\\.gz$", Thermo_Hela_files2))

        }
    }
    invisible(ans)
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
    .Defunct(new = "pxfile",
             package = "rpx",
             msg = paste0("Please use the rpx package to download data ",
                          "from the ProteomeXchaneg repository."))
}


##' Unless already present, downloads the PXD000001 mzTab file
##' in the \code{destdir} directory. The resulting file is named
##' \code{F063721.dat-mztab.txt}.
##'
##' @title Download the PXD000001 mzTab file
##' @param destdir A \code{character} with the destination folder.
##' @return Invisibly returns the name of the downloaded file.
##' @author Laurent Gatto
getPXD000001mzTab <- function(destdir = ".") {
    .Defunct(new = "pxfile",
             package = "rpx",
             msg = paste0("Please use the rpx package to download data ",
                          "from the ProteomeXchaneg repository."))
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
    .Defunct(new = "pxfile",
             package = "rpx",
             msg = paste0("Please use the rpx package to download data ",
                          "from the ProteomeXchaneg repository."))
}



