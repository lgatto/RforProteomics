.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This is the 'RforProteomics' version ",
                               packageVersion("RforProteomics"), ".\n",
                              "Run 'RforProtemics()' to get started.\n", sep=""))
  addVigs2WinMenu("synapter")
}

