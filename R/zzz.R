.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This is the 'RforProteomics' version ",
                               packageVersion("RforProteomics"), ".\n",
                              "Run 'RforProtemics()' or visit 'http://lgatto.github.com/RforProteomics/' to get started.\n", sep=""))
  addVigs2WinMenu("synapter")
}

