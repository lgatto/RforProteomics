\name{shinyMA}
\alias{shinyMA}


\title{
  MA and expression plots in \code{shiny}
}

\description{
  Starts an interactive \code{shiny} application that displays an MA
  plot and an expression plot side by side. The user can select features
  of interest on the MAplot and the respective intensities are displayed
  on the expression plot on the right.

  The data has been prepared using the mzTab file from the
  ProteomeXchange spiked-in data PXD000001 (see \code{\link{qnt}} for
  details). Sample 1, 2, 6 and 3, 4, 5 have been arbitrarily chosen to
  define two groups.  }

\usage{
shinyMA()
}


\value{
  Used for its side effects of starting a \code{shiny} application.
}

\references{
  The application is an adaptation of Michael Love's \code{shinyMA} app
  available on https://github.com/mikelove/shinyMA.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
if (interactive())
    shinyMA()
}

