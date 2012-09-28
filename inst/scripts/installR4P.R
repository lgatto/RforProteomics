installR4P <- function(v = "0.1.0") {
  suppressPackageStartupMessages(library("devtools"))
  message("[R4P] Installing dependencies...")
  deps <- c("R.utils", "Biobase",
            "mzR", "MSnbase", "xcms", "msdata", "isobar", 
            "MALDIquant", "readBrukerFlexData", 
            "synapter", "synapterdata",
            "IPPD", "Rdisop", "OrgMassSpecR", "BRAIN",
            "rols", "hpar", "GO.db", "org.Hs.eg.db", "biomaRt", 
            "RColorBrewer", "ggplot2", "reshape2",
            "knitr")
  instld <- rownames(installed.packages())
  deps <- deps[!deps %in% instld]
  if (length(deps) > 0) {
    if (!require("BiocInstaller")) {
      source("http://www.bioconductor.org/biocLite.R")
      useDevel(TRUE)
    }
    biocLite(deps, suppressUpdates = TRUE)
  }

  message("[R4P] Installing RforProteomics...")
  os <- .Platform$OS.type
  if (os == "unix") {
    ext <- ".tar.gz"
  } else {
    ext <- ".zip"
  }
  r4p <- paste0("RforProteomics_", v, ext)  
  url <- "https://github.com/downloads/lgatto/RforProteomics/"
  pkg <- paste0(url, r4p)

  if (os == "unix") {
    install_url(pkg)
  } else {
    if (!setInternet2(use = NA))
      setInternet2(use = TRUE)
    tdir <- tempdir()
    dest <- file.path(tdir, r4p)
    download.file(pkg, dest)
    install.packages(dest, type = "win.binary")
  }
  
  message("You can now type library('RforProteomics') to load the package.")
}

installR4P()


