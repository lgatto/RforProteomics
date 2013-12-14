\name{getThermoHelaPRTC}
\alias{getThermoHelaPRTC}
\title{Dowload Thermo Hela PRTC data}
\usage{
  getThermoHelaPRTC(src, destdir = ".")
}
\arguments{
  \item{src}{The name of the file to be downloaded. If
  missing, a vector of possible filenames is returned. If
  \code{"all"}, all files are downloaded. Alternatively, a
  pattern can be used to \code{grep} the files from the
  output \code{getThermoHelaPRTC()} the files to be
  downloaded.}

  \item{destdir}{Destination directory. Default is
  \code{"."}.}
}
\value{
  Invisibly return the path of the downloaded files.
}
\description{
  Downloads on of multiple Thermo Hela/PRTC data files.
}
\examples{
getThermoHelaPRTC()
getThermoHelaPRTC("design")
\dontrun{
getThermoHelaPRTC("all")
}
}
\author{
  Laurent Gatto
}
\seealso{
  \code{downloadData}
}

