installR4P <- function(v = "0.1.0") {
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
  
  message("[R4p] Downloading RforProteomics...")
  os <- .Platform$OS.type
  if (os == "unix") {
    ext <- ".tar.gz"
    type <- "source"
  } else {
    ext <- ".zip"
    type <- "win.binary"
  }
  r4p <- paste0("RforProteomics_", v, ext)
  
  url <- "https://github.com/downloads/lgatto/RforProteomics/"
  tdir <- tempdir()
  dest <- file.path(tdir, r4p)
  src <- paste0(url, r4p)

  download.file(url = src, destfile = dest)
  
  
  message("[R4P] Installing RforProteomics...")
  install.packages(dest, repos = NULL, type = type)
  message("You can now type library('RforProteomics') to load the package.")
}

installR4P()


