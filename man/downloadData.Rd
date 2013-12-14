\name{downloadData}
\alias{downloadData}
\title{Download a file}
\usage{
  downloadData(destdir = ".", src, unpack = TRUE, ...)
}
\arguments{
  \item{destdir}{The destination directory. Default is
  \code{"."}.}

  \item{src}{The url of the file to download.}

  \item{unpack}{Should \code{src} be uncompressed? Default
  is \code{TRUE}.}

  \item{...}{Additional paramters passed to
  \code{\link{download.file}}.}
}
\value{
  Invisible returns the full path of the downloaded file.
}
\description{
  Unless already present, downloads \code{src} in the
  \code{destdir} directory.
}
\author{
  Laurent Gatto
}

