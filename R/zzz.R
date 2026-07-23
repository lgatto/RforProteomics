.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("\nThis is the 'RforProteomics' version ",
                               utils::packageVersion("RforProteomics"), ".\n\n",
                               "  The material in this package is now outdated.\n",
                               "  Please visit the  R for Mass Spectrometry initiative\n",
                               "  website (https://www.rformassspectrometry.org/) and\n",
                               "  book (https://rformassspectrometry.github.io/book/)\n",
                               "  for the latest mass spectrometry and proteomics infrastructure.\n\n"))
}
