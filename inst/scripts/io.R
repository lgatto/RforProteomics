library(RforProteomics)
library(mzR)
library(MSnbase)

mzxml <- getPXD000001mzXML()
mzxml
raw <- openMSfile(mzxml)
raw

