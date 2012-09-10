library(mzR)
library(MSnbase)
library(R.utils)

download.file("ftp://ftp.pride.ebi.ac.uk/2012/03/PXD000001/PRIDE_Exp_Complete_Ac_22134.xml.gz",
              destfile = "PRIDE_Exp_Complete_Ac_22134.xml.gz")


gunzip("PRIDE_Exp_Complete_Ac_22134.xml.gz", unzip = "unzip")

raw <- openMSfile("PRIDE_Exp_Complete_Ac_22134.xml")
