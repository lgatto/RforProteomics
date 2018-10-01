.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("\nThis is the 'RforProteomics' version ",
                               utils::packageVersion("RforProteomics"), ".\n\n",
                               "  To get started, visit\n",
                               "    http://lgatto.github.com/RforProteomics/\n\n",
                               "  or, in R, open package vignettes by typing\n",
                               "    RforProteomics() # R/Bioc for proteomics overview\n",
                               "    RProtVis()       # R/Bioc for proteomics visualisation\n\n",
                               "  For a full list of available documents:\n",
                               "    vignette(package='RforProteomics')\n\n"))
}
