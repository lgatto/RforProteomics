v <- c("2.13", "2.14", "3.0", "3.1")
names(v) <- v

lpp <- lapply(v, RforProteomics::proteomicsPackages)
lmsp <- lapply(v, RforProteomics::massSpectrometryPackages)
lmsdp <- lapply(v, RforProteomics::massSpectrometryDataPackages)

saveRDS(lpp, file = "../extdata/lpp.rds")
saveRDS(lmsp, file = "../extdata/lmsp.rds")
saveRDS(lmsdp, file = "../extdata/lmsdp.rds")
