.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This is the 'RforProteomics' version ",
                               utils::packageVersion("RforProteomics"), ".\n",
                               "Run 'RforProtemics()' in R or visit \n",
                               "'http://lgatto.github.com/RforProteomics/' to get started.\n\n",
                               sep=""))
  addVigs2WinMenu("synapter")
}

