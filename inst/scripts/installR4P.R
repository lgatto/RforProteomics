installR4P <- function(v = "0.2.2", deps = TRUE) {
  stop("Please use 'BiocManager::install(\"RforProteomics\")' to install the package.")
  if (deps) {
    deps <- c("R.utils", "Biobase",
              "mzR", "MSnbase", "xcms", "msdata", "isobar", 
              "MALDIquant", "readBrukerFlexData", 
              "synapter", "synapterdata",
              "IPPD", "Rdisop", "OrgMassSpecR", "BRAIN",
              "rols", "hpar", "GO.db", "org.Hs.eg.db", "biomaRt", 
              "RColorBrewer", "ggplot2", "reshape2",
              "knitr")
  } else {
    deps <- c("R.utils", "Biobase")
  }    
  instld <- rownames(installed.packages())
  deps <- deps[!deps %in% instld]
  if (length(deps) > 0) {
    message("[R4P] Installing dependencies...")        
    if (!require("BiocInstaller")) {
      if (!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager")
      useDevel(TRUE)
    }
    BiocManager::install(deps, suppressUpdates = TRUE)
  }

  message("[R4P] Installing RforProteomics...")
  os <- .Platform$OS.type
  if (os == "unix") {
    ext <- ".tar.gz"
  } else {
    ext <- ".zip"
  }
  r4p <- paste0("RforProteomics_", v, ext)  
  url <- "http://proteome.sysbiol.cam.ac.uk/lgatto/RforProteomics/"
  pkg <- paste0(url, r4p)

  tdir <- tempdir()
  dest <- file.path(tdir, r4p)

  download.file(pkg, dest)
  
  if (os == "unix") {
    install.packages(dest, type = "source")
  } else {
    install.packages(dest, type = "win.binary")
  }
  
  message("You can now type library('RforProteomics') to load the package.")
}

installR4P()


