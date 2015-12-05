## ----biocstyle, eval=TRUE, echo=FALSE, results="asis"-----------------------------------
BiocStyle::latex()

## ----env0, echo = FALSE, warning = FALSE, message=FALSE---------------------------------
suppressPackageStartupMessages(library("BiocInstaller"))
## suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("MSnbase"))
suppressPackageStartupMessages(library("mzID"))
suppressPackageStartupMessages(library("rpx"))
suppressPackageStartupMessages(library("isobar"))
suppressPackageStartupMessages(library("MALDIquant"))
suppressPackageStartupMessages(library("MALDIquantForeign"))
suppressPackageStartupMessages(library("IPPD"))
suppressPackageStartupMessages(library("rols"))
suppressPackageStartupMessages(library("hpar"))
suppressPackageStartupMessages(library("BRAIN"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GO.db"))
suppressPackageStartupMessages(library("Rdisop"))
suppressPackageStartupMessages(library("rTANDEM"))
suppressPackageStartupMessages(library("MSGFplus"))
suppressPackageStartupMessages(library("MSGFgui"))

## ----'setup', include = FALSE, cache = FALSE--------------
library("knitr")
opts_chunk$set(tidy.opts =
                   list(width.cutoff = 50,
                   tidy = FALSE),
               fig.align = 'center',
               stop_on_error = 0L)
options(width = 60)

## ----vignette1, echo = TRUE, eval = FALSE-----------------
#  ## list all the vignettes in the RforProteomics package
#  vignette(package = "RforProteomics")
#  ## Open the vignette called RforProteomics
#  vignette("RforProteomics", package = "RforProteomics")
#  ## or just
#  vignette("RforProteomics")

## ----installation, eval = FALSE---------------------------
#  ## only first time you install Bioconductor packages
#  source("http://www.bioconductor.org/biocLite.R")
#  ## else
#  library("BiocInstaller")
#  biocLite("RforProteomics")

## ----installation2, eval = FALSE--------------------------
#  biocLite("RforProteomics", dependencies = TRUE)

## ----loadR4Prot, warning=FALSE----------------------------
library("RforProteomics")

## ----stangle, eval=TRUE, tidy=FALSE, error=FALSE----------
## gets the vignette source
rnwfile <- system.file("doc/vigsrc/RforProteomics.Rnw",
                       package = "RforProteomics")
## produces the R file in the working directory
library("knitr")
purl(rnwfile, quiet = TRUE)

## ----env, message=FALSE-----------------------------------
library("RColorBrewer") ## Color palettes
library("ggplot2")  ## Convenient and nice plotting
library("reshape2") ## Flexibly reshape data

## ----mzr--------------------------------------------------
## load the required packages
library("mzR") ## the software package
library("msdata") ## the data package
## below, we extract the releavant example file
## from the local 'msdata' installation
filepath <- system.file("microtofq", package = "msdata")
file <- list.files(filepath, pattern="MM14.mzML",
                   full.names=TRUE, recursive = TRUE)
## creates a commection to the mzML file
mz <- openMSfile(file)
## demonstraction of data access
basename(fileName(mz))
isInitialized(mz)
runInfo(mz)
instrumentInfo(mz)
## once finished, it is good to explicitely
## close the connection
close(mz)

## ----mzrid------------------------------------------------
file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
mzid <- openIDfile(file)
mzid

## ----mzrid2-----------------------------------------------
softwareInfo(mzid)
enzymes(mzid)
names(psms(mzid))
head(psms(mzid))[, 1:13]

## ----mzid1------------------------------------------------
library("mzID")
id <- mzID("http://psi-pi.googlecode.com/svn/trunk/examples/1_1examples/55merge_tandem.mzid")
id

## ----mzid2------------------------------------------------
ids <- mzID(c(
    "http://psi-pi.googlecode.com/svn/trunk/examples/1_1examples/55merge_tandem.mzid",
    "http://psi-pi.googlecode.com/svn/trunk/examples/1_1examples/55merge_omssa.mzid"))
ids

## ----flatid-----------------------------------------------
fid <- flatten(id)
names(fid)
dim(fid)

## ----msnexp-----------------------------------------------
library("MSnbase")
## uses a simple dummy test included in the package
mzXML <- dir(system.file(package="MSnbase",dir="extdata"),
             full.name=TRUE,
             pattern="mzXML$")
basename(mzXML)
## reads the raw data into and MSnExp instance
raw <- readMSData(mzXML, verbose = FALSE)
raw
## Extract a single spectrum
raw[[3]]

## ----msnexpPlot, dev='pdf', echo=TRUE, fig.width=5, fig.height=5, fig.keep='all', out.width='.48\\linewidth'----
plot(raw, full=TRUE)
plot(raw[[3]], full=TRUE, reporters=iTRAQ4)

## ----downloadmztab, tidy = FALSE--------------------------
## Experiment information
library("rpx")
px1 <- PXDataset("PXD000001")
px1
pxfiles(px1)
## Downloading the mzTab data
mztab <- pxget(px1, "PXD000001_mztab.txt")
mztab

## ----mztab, tidy = FALSE----------------------------------
## Load mzTab peptide data
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")
sampleNames(qnt) <- reporterNames(TMT6)
head(exprs(qnt))
## remove missing values
qnt <- filterNA(qnt)
processingData(qnt)

## combine into proteins
## - using the 'accession' feature meta data
## - sum the peptide intensities
protqnt <- combineFeatures(qnt,
                           groupBy = fData(qnt)$accession,
                           fun = sum)

## ----matplot, dev='pdf', echo=TRUE, fig.width=5.5, fig.height=5.5, fig.keep='last', out.width='.6\\linewidth', tidy=FALSE----
cls <- brewer.pal(5, "Set1")
matplot(t(tail(exprs(protqnt), n = 5)), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(protqnt), n=5),
       lty = 1, bty = "n", cex = .8, col = cls)

## ----mztab2-----------------------------------------------
qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")

acc <- c("P00489", "P00924",
         "P02769", "P62894",
         "ECA")

idx <- sapply(acc, grep, fData(qnt)$accession)
idx2 <- sapply(idx, head, 3)
small <- qntS[unlist(idx2), ]

idx3 <- sapply(idx, head, 10)
medium <- qntV[unlist(idx3), ]

m <- exprs(medium)
colnames(m) <- c("126", "127", "128",
                 "129", "130", "131")
rownames(m) <- fData(medium)$accession
rownames(m)[grep("CYC", rownames(m))] <- "CYT"
rownames(m)[grep("ENO", rownames(m))] <- "ENO"
rownames(m)[grep("ALB", rownames(m))] <- "BSA"
rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
rownames(m)[grep("ECA", rownames(m))] <- "Background"

cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),
         "grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)

## ----heatmap, dev='pdf', echo=TRUE, fig.width=5.5, fig.height=5.5, fig.keep='last', out.width='.6\\linewidth'----
heatmap(m, col = wbcol, RowSideColors=cls[rownames(m)])

## ----spikes, dev='pdf', echo=TRUE, fig.width=9, fig.height=5, fig.keep='all', out.width='1\\linewidth', tidy=FALSE----
dfr <- data.frame(exprs(small),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)
colnames(dfr) <- c("126", "127", "128", "129", "130", "131",
                   "Protein", "Feature")
dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr2 <- melt(dfr)
ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity")

## ----qntdf------------------------------------------------
d <- data.frame(Signal = rowSums(exprs(qntms)[, 1:6]),
                Incomplete = exprs(qntms)[, 7])
d <- log(d)
cls <- rep("#00000050", nrow(qnt))
pch <- rep(1, nrow(qnt))
cls[grep("P02769", fData(qnt)$accession)] <- "gold4" ## BSA
cls[grep("P00924", fData(qnt)$accession)] <- "dodgerblue" ## ENO
cls[grep("P62894", fData(qnt)$accession)] <- "springgreen4" ## CYT
cls[grep("P00489", fData(qnt)$accession)] <- "darkorchid2" ## PHO
pch[grep("P02769", fData(qnt)$accession)] <- 19
pch[grep("P00924", fData(qnt)$accession)] <- 19
pch[grep("P62894", fData(qnt)$accession)] <- 19
pch[grep("P00489", fData(qnt)$accession)] <- 19

## ----plotmzdelta, dev='pdf', echo=TRUE, fig.width=8, fig.height=5, fig.keep='last', out.width='0.9\\linewidth', warning = FALSE----
mzp

## ----incompl, dev='pdf', echo=TRUE, fig.width=5.5, fig.height=5.5, fig.keep='last', out.width='.6\\linewidth', tidy=FALSE----
plot(Signal ~ Incomplete, data = d,
     xlab = expression(Incomplete~dissociation),
     ylab = expression(Sum~of~reporters~intensities),
     pch = 19,
     col = "#4582B380")
grid()
abline(0, 1, lty = "dotted")
abline(lm(Signal ~ Incomplete, data = d), col = "darkblue")

## ----maplot, dev='pdf', echo=TRUE, fig.width=5.5, fig.height=5.5, fig.keep='last', out.width='.6\\linewidth'----
MAplot(qnt[, c(4, 2)], cex = .9, col = cls, pch = pch, show.statistics = FALSE)

## ----mqload, tidy=FALSE-----------------------------------
## load packages
library("MALDIquant")
library("MALDIquantForeign")
## getting test data
datapath <-
  file.path(system.file("Examples",
                        package = "readBrukerFlexData"),
            "2010_05_19_Gibb_C8_A1")
dir(datapath)
sA1 <- importBrukerFlex(datapath, verbose=FALSE)
# in the following we use only the first spectrum
s <- sA1[[1]]

summary(mass(s))
summary(intensity(s))
head(as.matrix(s))

## ----mqplot, dev='pdf', echo=TRUE, fig.width=7, fig.height=5.5, fig.keep='last', out.width='.6\\linewidth'----
plot(s)

## ----mqpreproc--------------------------------------------
## sqrt transform (for variance stabilization)
s2 <- transformIntensity(s, method="sqrt")
s2

## smoothing - 5 point moving average
s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
s3

## baseline subtraction
s4 <- removeBaseline(s3, method="SNIP")
s4

## ----mqred------------------------------------------------
## peak picking
p <- detectPeaks(s4)
length(p) # 181
peak.data <- as.matrix(p) # extract peak information

## ----mqplot2, dev='pdf', echo=TRUE, fig.width=10, fig.height=5, fig.keep='high', out.width='1\\linewidth'----
par(mfrow=c(2,3))
xl <- range(mass(s))
# use same xlim on all plots for better comparison
plot(s, sub="", main="1: raw", xlim=xl)
plot(s2, sub="", main="2: variance stabilisation", xlim=xl)
plot(s3, sub="", main="3: smoothing", xlim=xl)
plot(s4, sub="", main="4: base line correction", xlim=xl)
plot(s4, sub="", main="5: peak detection", xlim=xl)
points(p)
top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
labelPeaks(p, index=top20, underline=TRUE)
plot(p, sub="", main="6: peak plot", xlim=xl)
labelPeaks(p, index=top20, underline=TRUE)

## ----isotopes, tidy = FALSE, warning = FALSE--------------
library(IPPD)
library(BRAIN)
atoms <- getAtomsFromSeq("SIVPSGASTGVHEALEMR")
unlist(atoms)

library(Rdisop)
pepmol <- getMolecule(paste0(names(atoms),
                             unlist(atoms),
                             collapse = ""))
pepmol

##
library(OrgMassSpecR)
data(itraqdata)

simplottest <-
  itraqdata[featureNames(itraqdata) %in% paste0("X", 46:47)]
sim <- SpectrumSimilarity(as(simplottest[[1]], "data.frame"),
                          as(simplottest[[2]], "data.frame"),
                          top.lab = "itraqdata[['X46']]",
                          bottom.lab = "itraqdata[['X47']]",
                          b = 25)
title(main = paste("Spectrum similarity", round(sim, 3)))

MonoisotopicMass(formula = list(C = 2, O = 1, H=6))
molecule <- getMolecule("C2H5OH")
molecule$exactmass
## x11()
## plot(t(.pepmol$isotopes[[1]]), type = "h")

## x <- IsotopicDistribution(formula = list(C = 2, O = 1, H=6))
## t(molecule$isotopes[[1]])
## par(mfrow = c(2,1))
## plot(t(molecule$isotopes[[1]]), type = "h")
## plot(x[, c(1,3)], type = "h")

## data(myo500)
## masses <- c(147.053, 148.056)
## intensities <- c(93, 5.8)
## molecules <- decomposeIsotopes(masses, intensities)

## experimental eno peptides
exppep <-
  as.character(fData(qnt[grep("ENO", fData(qnt)[, 2]), ])[, 1]) ## 13
minlength <- min(nchar(exppep))


if (!file.exists("P00924.fasta"))
    eno <- download.file("http://www.uniprot.org/uniprot/P00924.fasta",
                         destfile = "P00924.fasta")
eno <- paste(readLines("P00924.fasta")[-1], collapse = "")
enopep <- Digest(eno, missed = 1)
nrow(enopep) ## 103
sum(nchar(enopep$peptide) >= minlength) ## 68
pepcnt <- enopep[enopep[, 1] %in% exppep, ]
nrow(pepcnt) ## 13

## ----cleaver, tidy = FALSE--------------------------------
library(cleaver)
cleave("LAAGKVEDSD", enzym = "trypsin")

## ----cleaver_missing, tidy = FALSE------------------------
## miss one cleavage position
cleave("LAAGKVEDSD", enzym = "trypsin", missedCleavages = 1)

## miss zero or one cleavage positions
cleave("LAAGKVEDSD", enzym = "trypsin", missedCleavages = 0:1)

## ----ibplot, dev='pdf', echo=TRUE, fig.width=5, fig.height=5, fig.keep='last'----
maplot(x,
       noise.model = c(nm.background, nm.spks, nm),
       channel1="127", channel2="129",
       pch = 19, col = cls2,
       main = "Spectra MA plot")
abline(h = 1, lty = "dashed", col = "grey")
legend("topright",
       c("BSA", "ENO", "CYT", "PHO"),
       pch = 19, col = c("gold4", "dodgerblue",
                   "springgreen4", "darkorchid2"),
       bty = "n", cex = .7)

## ----synapter, eval=FALSE---------------------------------
#  ## open the synapter vignette
#  library("synapter")
#  synapterGuide()

## ----rtandem1, tidy = FALSE-------------------------------
library(rTANDEM)
taxonomy <- rTTaxo(taxon = "yeast",
                   format = "peptide",
                   URL = system.file(
                     "extdata/fasta/scd.fasta.pro",
                     package="rTANDEM"))
param <- rTParam()
param <- setParamValue(param,
                       'protein', 'taxon',
                       value="yeast")
param <- setParamValue(param, 'list path',
                       'taxonomy information', taxonomy)
param <- setParamValue(param,
                       'list path', 'default parameters',
                       value = system.file(
                         "extdata/default_input.xml",
                         package="rTANDEM"))
param <- setParamValue(param, 'spectrum', 'path',
                       value = system.file(
                         "extdata/test_spectra.mgf",
                         package="rTANDEM"))
param <- setParamValue(param, 'output', 'xsl path',
                       value = system.file(
                         "extdata/tandem-input-style.xsl",
                         package="rTANDEM"))
param <- setParamValue(param, 'output', 'path',
                       value = paste(getwd(),
                         "output.xml", sep="/"))

## ----rtandem2---------------------------------------------
resultPath <- tandem(param)
basename(resultPath)

## ----rtandem3---------------------------------------------
res <- GetResultsFromXML(resultPath)
## the inferred proteins
proteins <- GetProteins(res,
                        log.expect = -1.3,
                        min.peptides = 2)
proteins[, -(4:5), with = FALSE]
## the identified peptides for YFR053C
peptides <- GetPeptides(protein.uid = 1811,
                        results = res,
                        expect = 0.05)
peptides[, c(1:4, 9, 10:16), with = FALSE]

## ----msgf1------------------------------------------------
library("MSGFplus")
## Create a parameter object with a set of parameters
param <- msgfPar(database = system.file('extdata',
                                        'milk-proteins.fasta',
                                        package='MSGFplus'),
                 tolerance = '10 ppm',
                 enzyme = 'Trypsin')

## Add parameters after creation
instrument(param) <- 'QExactive'
tda(param) <- TRUE
ntt(param) <- 2

## Add expected modifications
mods(param)[[1]] <- msgfParModification('Carbamidomethyl',
                                        composition = 'C2H3N1O1',
                                        residues = 'C',
                                        type = 'fix',
                                        position = 'any')

mods(param)[[2]] <- msgfParModification(name = 'Oxidation',
                                        mass = 15.994915,
                                        residues = 'M',
                                        type = 'opt',
                                        position = 'any')

nMod(param) <- 2    # Number of allowed modifications per peptide

## Get a summary of your parameters
show(param)

## ----msgfplus2, eval=FALSE--------------------------------
#  result <- runMSGF(param, 'path/to/a/rawfile.mzML')

## ----MSnIDstart-------------------------------------------
library("MSnID")
msnid <- MSnID(".")

## ----MSnIDdataImport1-------------------------------------
PSMresults <- read.delim(system.file("extdata", "human_brain.txt",
                                     package="MSnID"),
                         stringsAsFactors=FALSE)
psms(msnid) <- PSMresults
show(msnid)

## ----MSnIDsequence----------------------------------------
msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
prop.table(table(msnid$numIrregCleavages))

## ----MSnIDmissedCleavagesPlot, dev='pdf', fig.width=8, fig.height=4, out.width='.7\\textwidth'----
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) +
    geom_bar(stat='identity', position='dodge') +
    ggtitle("Number of Missed Cleavages")

## ----MSnIDfilteringCriteria-------------------------------
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

## ----MSnIDsettingFilter-----------------------------------
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
filtObj$msmsScore <- list(comparison=">", threshold=8.0)
show(filtObj)

## ----MSnIDfilterAssessment--------------------------------
evaluate_filter(msnid, filtObj, level="PSM")
evaluate_filter(msnid, filtObj, level="peptide")
evaluate_filter(msnid, filtObj, level="accession")

## ----MSnIDfilterOptimization1-----------------------------
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)
show(filtObj.grid)

## ----MSnIDfilterOptimization2-----------------------------
filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01,
                                method="Nelder-Mead", level="peptide",
                                n.iter=500)
show(filtObj.nm)

## ----MSnIDfilterAssessment2-------------------------------
evaluate_filter(msnid, filtObj, level="peptide")
evaluate_filter(msnid, filtObj.grid, level="peptide")
evaluate_filter(msnid, filtObj.nm, level="peptide")

## ----MSnIDapplyFilter-------------------------------------
msnid <- apply_filter(msnid, filtObj.nm)
show(msnid)

## ----MSnIDremovingDecoyAndContaminants--------------------
msnid <- apply_filter(msnid, "isDecoy == FALSE")
show(msnid)
msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
show(msnid)

## ----MSnIDgetDataOut1-------------------------------------
psm.df <- psms(msnid)
psm.dt <- as(msnid, "data.table")

## ----MSnIDgetDataOut2-------------------------------------
peps <- peptides(msnid)
head(peps)
prots <- accessions(msnid)
head(prots)
prots <- proteins(msnid) # may be more intuitive then accessions
head(prots)

## ----MSnIDconvertingToMSnSet------------------------------
msnset <- as(msnid, "MSnSet")
library("MSnbase")
head(fData(msnset))
head(exprs(msnset))

## ----MSnIDsummingPeptidesToProteins-----------------------
msnset <- combineFeatures(msnset,
                            fData(msnset)$accession,
                            redundancy.handler="unique",
                            fun="sum",
                            cv=FALSE)
head(fData(msnset))
head(exprs(msnset))

## ----eval=TRUE, echo=FALSE, results='hide'----------------
unlink(".Rcache", recursive=TRUE)

## ----idex1, message=FALSE---------------------------------
library("rpx")
id <- "PXD002161"
px <- PXDataset(id)
try(setInternet2(FALSE),silent=TRUE)
library("jsonlite")
addr <- "http://www.ebi.ac.uk:80/pride/ws/archive/%s/list/project/%s"
files <- fromJSON(sprintf(addr, "file", id))$list
assays <- fromJSON(sprintf(addr, "assay", id))$list

files <- subset(files, fileType == 'PEAK',
                select = c("assayAccession","fileName"))
assays <- assays[,c("assayAccession",
                    "experimentalFactor",
                    "proteinCount",
                    "peptideCount",
                    "uniquePeptideCount",
                    "identifiedSpectrumCount",
                    "totalSpectrumCount")]

pttrn <- "Age: young, Diet; fully fed"
assays <- subset(assays, grepl(pttrn, experimentalFactor, fixed=T))
# encode phenotype and sample names
group <- sub(".*Name: Y-(.+?)-FF\\.(\\d)", "\\1", assays$experimentalFactor)
splnm <- sub(".*Name: Y-(.+?)-FF\\.(\\d)", "\\1_\\2", assays$experimentalFactor)

assays <- with(assays, {data.frame(assayAccession,
                                   phenotype=sub(".*Name: Y-(.+?)-FF\\.(\\d)",
                                                 "\\1", experimentalFactor),
                                   sampleName=sub(".*Name: Y-(.+?)-FF\\.(\\d)", 
                                                  "\\1_\\2", experimentalFactor),
                                   stringsAsFactors=F)})
files <- subset(files, assayAccession %in% assays$assayAccession)

## do we already have all files?
allfiles <- length(dir(pattern = "c_elegans_.*mzML")) == 10

if (!allfiles) 
    pxget(px, files$fileName) # fetch them from PX

files$datasetName <- sub('.mzML.gz','', files$fileName, fixed=TRUE)
meta <- merge(files[,c("assayAccession","datasetName")], assays)
rownames(meta) <- meta$datasetName
meta <- meta[order(meta$sampleName),]
rownames(meta) <- NULL

if (!allfiles) {
    library("R.utils")
    sapply(list.files(pattern = "mzML.gz"), gunzip)
}

## ----idex2, eval=FALSE------------------------------------
#  library("Biostrings")
#  fasta_location <-  "ftp://ftp.wormbase.org/pub/wormbase/releases/WS250/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS250.protein.fa.gz"
#  fwd.seqs <- readAAStringSet(fasta_location, format="fasta",
#                              nrec=-1L, skip=0L, use.names=TRUE)
#  rev.seqs <- reverse(fwd.seqs)
#  names(rev.seqs) <- paste("XXX", names(rev.seqs), sep='_')
#  fwd.rev.seqs <- append( fwd.seqs, rev.seqs)
#  writeXStringSet(x=fwd.rev.seqs, filepath="c_elegans_fwd_rev.fasta", format="fasta")

## ----idex3------------------------------------------------
library("rTANDEM")
param <- setParamOrbitrap()
taxonomy <- rTTaxo(taxon="celegans",
                   format="peptide",
                   URL= "c_elegans_fwd_rev.fasta")
param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
param <- setParamValue(param, 'protein', 'taxon', value='celegans')

def.input.path <- system.file("extdata/default_input.xml", package="rTANDEM")
param <- setParamValue(param, 'list path', 'default parameters',
                       value=def.input.path)

param <- setParamValue(param, "output", "message", "r-for-proteomics ")
param <- setParamValue(param, "refine", value="no")

## ----expar------------------------------------------------
library("parallel")
param <- setParamValue(param, "spectrum", "threads", detectCores())

## ----rtout------------------------------------------------
output.files <- lapply(sub("\\.gz","",files$fileName),
                       function(x){
                           param <- setParamValue(param, 'spectrum', 'path', value=x)
                           output.file <- tandem(param)})

## ----xttodf, warning=FALSE, message=FALSE-----------------
library("dplyr")
.P <- MSnID:::.PROTON_MASS # 1.007276

## parse to comply with MSnID
xtandem_to_df <- function(output.file){
    res <- GetResultsFromXML(output.file)
    proteins <- GetProteins(res)
    peptides <- GetPeptides(protein.uid = proteins$uid, results = res)
    peptides <- select(as.data.frame(res@peptides), 
                       prot.uid, spectrum.mh, mh, expect.value, 
                       tandem.score, start.position, end.position,
                       spectrumID = spectrum.id,
                       chargeState = spectrum.z,
                       pepSeq = sequence)
    peptides <- mutate(peptides, 
                       calculatedMassToCharge = (mh + .P*(chargeState - 1))/chargeState,
                       experimentalMassToCharge = (spectrum.mh + .P*(chargeState - 1))/chargeState,
                       spectrum.mh = NULL,
                       mh = NULL)
    proteins <- select(as.data.frame(res@proteins),
                       prot.uid=uid, sequence, description)
    proteins <- mutate(proteins, 
                       accession = sub("(\\S+)\t.*", "\\1", description),
                       isDecoy = grepl("XXX_", accession))
    x <- inner_join(peptides, proteins)
    pre <- with(x, substr(sequence, start.position-1, start.position-1))
    pre <- ifelse(pre == "", "-", pre)
    post <- with(x, substr(sequence, end.position+1, end.position+1))
    post <- ifelse(post == "", "-", post)
    x <- mutate(x, 
                peptide = paste(pre,x$pepSeq, post, sep='.'),
                sequence = NULL)
    return(x)
}

out <- lapply(output.files, xtandem_to_df)
for(i in seq_along(out))
    out[[i]]$spectrumFile <- sub("\\.mzML\\.gz","",files$fileName)[i]
out <- Reduce(rbind, out)

## ----msnidm-----------------------------------------------
library(MSnID)
m <- MSnID()
psms(m) <- out
show(m)

## ----msnisfig, fig.width=6, fig.height=4------------------
library("ggplot2")
dM <- with(psms(m),
             round((experimentalMassToCharge-calculatedMassToCharge)*chargeState))
x <- as.data.frame(table(data.frame(dM, isDecoy=m$isDecoy)))
x <- as.data.frame(table(dM, isDecoy=m$isDecoy))
ggplot(x, aes(x=dM, fill=isDecoy, y=Freq)) +
      geom_bar(stat="identity", width=0.25)

## ----correctpeaksel, fig.with=6, fig.height=6-------------
m <- correct_peak_selection(m)
dM <- with(psms(m),
             round((experimentalMassToCharge-calculatedMassToCharge)*chargeState))
x <- as.data.frame(table(data.frame(dM, isDecoy=m$isDecoy)))
x <- as.data.frame(table(dM, isDecoy=m$isDecoy))
ggplot(x, aes(x=dM, fill=isDecoy, y=Freq)) +
      geom_bar(stat="identity", width=0.25)

## ----parentionmass, fig.width=8, fig.height=6-------------
mme <- data.frame(ppm=mass_measurement_error(m), isDecoy=m$isDecoy)
ggplot(mme, aes(x=ppm, fill=isDecoy)) +
      geom_histogram(binwidth=1) +
      xlim(-20,+20) +
      facet_wrap( ~ isDecoy) + 
   	scale_y_sqrt()

## ----msnidfilt--------------------------------------------
m$score <- -log10(m$expect.value) # higher the better
m$mme <- abs(mass_measurement_error(m)) # lower the better
fltr <- MSnIDFilter(m)

## ---------------------------------------------------------
fltr$score <- list(comparison=">", threshold=0)
fltr$mme <- list(comparison="<", threshold=0)
fltr.grid <- optimize_filter(fltr, m, fdr.max=0.01,
method="Grid", level="peptide",
n.iter=1000)
show(fltr.grid)

## ---------------------------------------------------------
as.numeric(fltr.grid)
evaluate_filter(m, fltr.grid)

## ----applsann---------------------------------------------
m <- apply_filter(m, fltr.sann)

## ----nbidpeps, fig.width=6, fig.height=6------------------
pp <- unique(psms(m)[,c("peptide","accession")])
pep_per_prot <- as.data.frame(table(table(pp$accession)))
pep_per_prot$peptides <- with(pep_per_prot, 
                              ifelse(as.numeric(Var1) > 5, "6+", as.numeric(Var1)))
ggplot(pep_per_prot, 
         aes(x=factor(1), fill=peptides, y=Freq)) + 
      geom_bar(width=1, stat = "identity", position = "stack") + 
      coord_polar(theta = "y") + 
      ylab('') + xlab('') + 
      scale_fill_discrete(name="number of\npeptides per\nprotein")

## ----applflitbigex----------------------------------------
nms <- names(which(table(pp$accession) > 1))
m <- apply_filter(m, "accession %in% nms")
show(m)

## ---------------------------------------------------------
mset <- as(m, "MSnSet")
show(mset)

## ---------------------------------------------------------
library("MSnbase")
mset <- combineFeatures(mset,
fData(mset)$accession,
redundancy.handler="unique",
fun="sum",
cv=FALSE)
show(mset)

## ---------------------------------------------------------
mset <- mset[rowSums(exprs(mset) > 0) >= 6,] 

## ---------------------------------------------------------
pData(mset) <- meta[match(sampleNames(mset), meta$datasetName),]
mset$phenotype <- as.factor(mset$phenotype)
sampleNames(mset) <- mset$sampleName
pData(mset)

## ----scpca, fig.width=6, fig.height=6---------------------
library("msmsEDA")
counts.pca(mset,
             facs = pData(mset)[,'phenotype',drop=FALSE],
             snms = sampleNames(mset))

## ----phist, fig.height=6, fig.width=7---------------------
library("msmsTests")
  alt.f <- "y ~ phenotype"
    null.f <- "y ~ 1"
      div <- colSums(exprs(mset)) # normalization factor
 	res <- msms.glm.qlll(mset, "y ~ phenotype", "y ~ 1", div=div)
 	hist(res$p.value, 100)

## ----scvolc, fig.width=6, fig.height=6--------------------
lst <- test.results(res, mset, pData(mset)$phenotype,
                      "ctrl","daf2", div, alpha=0.05,
                      minSpC=0, minLFC=1, method="BH")
res.volcanoplot(lst$tres, min.LFC=1, max.pval=0.05, ylbls=NULL, maxy=4)

## ----hmap, fig.width=7, fig.height=7----------------------
library("Heatplus")
regulated <- subset(lst$tres, adjp < 0.05 & abs(LogFC) > 1)
## order MSnSet object the daf-16 status
mset <- mset[,order(pData(mset)$phenotype)]
## matrix with regulated proteins
selected.data <- exprs(mset[rownames(regulated),])
## scaling counts from 0 to 1
selected.data <- sweep(selected.data, 1, apply(selected.data, 1, min), '-')
selected.data <- sweep(selected.data, 1, apply(selected.data, 1, max), '/')
heatmap_plus(selected.data,
               scale='none',
               col=colorRampPalette(c("snow","steelblue"))(10))

## ----annot1, cache=FALSE----------------------------------
id <- "ENSG00000105323"
library("hpar")
getHpa(id, "SubcellularLoc")

## ----annot2, cache=FALSE, warning=FALSE-------------------
library("org.Hs.eg.db")
library("GO.db")
ans <- AnnotationDbi::select(org.Hs.eg.db,
keys = id,
columns = c("ENSEMBL", "GO", "ONTOLOGY"),
keytype = "ENSEMBL")
ans <- ans[ans$ONTOLOGY == "CC", ]
ans
sapply(as.list(GOTERM[ans$GO]), slot, "Term")

## ----xtabpkg, echo=FALSE, message=FALSE-------------------
library("xtable")

## ----protstab, echo=FALSE, results='asis'-----------------
t1 <- xtable(pkTab[["Proteomics"]],
             caption = "Packages available under the \\texttt{Proteomics} \\textit{biocViews} category.",
             table.placement = "thb",
             label = "tab:prot")
align(t1) <- c("l", "l", "p{10cm}", "l")
print(t1,
      include.rownames = FALSE,
      floating = FALSE,
      tabular.environment = 'longtable',
      size = "scriptsize")

## ----mstab, echo=FALSE, results='asis'--------------------
t2 <- xtable(pkTab[["MassSpectrometry"]],
             caption = "Packages available under the \\texttt{MassSpectrometry} \\textit{biocViews} category.",
             table.placement = "thb",
             label = "tab:ms")
align(t2) <- c("l", "l", "p{10cm}", "l")
print(t2,
      include.rownames = FALSE,
      floating = FALSE,
      tabular.environment = 'longtable',
      size = "scriptsize")

## ----msdatatab, echo=FALSE, results='asis'----------------
t3 <- xtable(pkTab[["MassSpectrometryData"]],
             caption = "Experimental Packages available under the \\texttt{MassSpectrometryData} \\textit{biocViews} category.",
             table.placement = "thb",
             label = "tab:msdata")
align(t3) <- c("l", "l", "p{10cm}", "l")
print(t3,
      include.rownames = FALSE,
      floating = FALSE,
      tabular.environment = 'longtable',
      usize = "scriptsize")

## ----pkgs, eval=FALSE-------------------------------------
#  pp <- proteomicsPackages()
#  display(pp)

## ----cntCRAN, echo=FALSE----------------------------------
X <- readLines("http://cran.r-project.org/web/views/ChemPhys.html")
x2 <- grep("Related links:", X)
x1 <- grep("CRAN packages:", X)
np <- length(grep("../packages/", X[x1:x2]))

## ----sessioninfo, results='asis', echo=FALSE, echo=FALSE----
toLatex(sessionInfo(), locale = FALSE)

