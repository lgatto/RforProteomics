## ----style, echo = FALSE, results = 'asis', message=FALSE----------------
BiocStyle::markdown()

## ----env, message=FALSE, echo=FALSE, warning=FALSE-----------------------
library("RforProteomics")
library("BiocInstaller")
library("mzR")
library("MSnbase")
library("knitr")
library("rpx")
library("xtable")
library("RColorBrewer")
library("MALDIquant")
library("MALDIquantForeign")
library("pRoloc")
library("pRolocdata")
library("msmsTests")
library("msmsEDA")

## ----packs, cache=FALSE, warning=FALSE, echo=FALSE-----------------------
library("RforProteomics")
pp <- proteomicsPackages()
msp <- massSpectrometryPackages()

## ----pptab, echo=FALSE, results='asis'-----------------------------------
kable(pp, format = "html")

## ----msptab, echo=FALSE, results='asis'----------------------------------
kable(msp, format = "html")

## ----anscombe, echo = FALSE, results='asis'------------------------------
kable(anscombe, format = "html")

## ----anscombetab---------------------------------------------------------

tab <- matrix(NA, 5, 4)
colnames(tab) <- 1:4
rownames(tab) <- c("var(x)", "mean(x)",
                   "var(y)", "mean(y)",
                   "cor(x,y)")

for (i in 1:4)
    tab[, i] <- c(var(anscombe[, i]),
                  mean(anscombe[, i]),
                  var(anscombe[, i+4]),
                  mean(anscombe[, i+4]),
                  cor(anscombe[, i], anscombe[, i+4]))


## ----anstabdisplay, echo=FALSE-------------------------------------------
kable(tab)

## ----anscombefig, dev='png'----------------------------------------------
ff <- y ~ x

mods <- setNames(as.list(1:4), paste0("lm", 1:4))

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
for (i in 1:4) {
    ff[2:3] <- lapply(paste0(c("y","x"), i), as.name)
    plot(ff, data = anscombe, pch = 19, xlim = c(3, 19), ylim = c(3, 13))
    mods[[i]] <- lm(ff, data = anscombe)
    abline(mods[[i]])
}

## ----anscomberesids, echo=FALSE, fig.cap="The 11 sets of residuals for Anscombe's four datasets."----
kable(sapply(mods, residuals))

## ----, makemadata, warning=FALSE, cache=FALSE----------------------------
library("rpx")
px1 <- PXDataset("PXD000001")
mztab <- pxget(px1, "PXD000001_mztab.txt")

library("MSnbase")
qnt <- readMzTabData(mztab, what = "PEP")
sampleNames(qnt) <- reporterNames(TMT6)
qnt <- filterNA(qnt)
## may be combineFeatuers

spikes <- c("P02769", "P00924", "P62894", "P00489")
protclasses <- as.character(fData(qnt)$accession)
protclasses[!protclasses %in% spikes] <- "Background"


madata42 <- data.frame(A = rowMeans(log(exprs(qnt[, c(4, 2)]), 10)),
                       M = log(exprs(qnt)[, 4], 2) - log(exprs(qnt)[, 2], 2),
                       data = rep("4vs2", nrow(qnt)),
                       protein = fData(qnt)$accession,
                       class = protclasses)

madata62 <- data.frame(A = rowMeans(log(exprs(qnt[, c(6, 2)]), 10)),
                       M = log(exprs(qnt)[, 6], 2) - log(exprs(qnt)[, 2], 2),
                       data = rep("6vs2", nrow(qnt)),
                       protein = fData(qnt)$accession,
                       class = protclasses)


madata <- rbind(madata42, madata62)

## ----, mafig1------------------------------------------------------------
par(mfrow = c(1, 2))
plot(M ~ A, data = madata42, main = "4vs2",
     xlab = "A", ylab = "M", col = madata62$class)
plot(M ~ A, data = madata62, main = "6vs2",
     xlab = "A", ylab = "M", col = madata62$class)


## ----mafig2--------------------------------------------------------------
library("lattice")
latma <- xyplot(M ~ A | data, data = madata,
                groups = madata$class,
                auto.key = TRUE)
print(latma)


## ----mafig3--------------------------------------------------------------

library("ggplot2")
ggma <- ggplot(aes(x = A, y = M, colour = class), data = madata,
               colour = class) +
                   geom_point() +
                       facet_grid(. ~ data)
print(ggma)


## ----macols--------------------------------------------------------------

library("RColorBrewer")
bcols <- brewer.pal(4, "Set1")
cls <- c("Background" = "#12121230",
         "P02769" = bcols[1],
         "P00924" = bcols[2],
         "P62894" = bcols[3],
         "P00489" = bcols[4])


## ----macust--------------------------------------------------------------

ggma2 <- ggplot(aes(x = A, y = M, colour = class),
                data = madata) + geom_point(shape = 19) +
                    facet_grid(. ~ data) + scale_colour_manual(values = cls) +
                        guides(colour = guide_legend(override.aes = list(alpha = 1)))
print(ggma2)


## ----mafigmsnset---------------------------------------------------------
MAplot(qnt, cex = .8)

## ----shinyMA, eval=FALSE-------------------------------------------------
## shinyMA()

## ----msmsTestsData, cache=FALSE------------------------------------------
library("msmsEDA")
library("msmsTests")
data(msms.dataset)
## Pre-process expression matrix
e <- pp.msms.data(msms.dataset)
## Models and normalizing condition
null.f <- "y~batch"
alt.f <- "y~treat+batch"
div <- apply(exprs(e),2,sum)
## Test
res <- msms.glm.qlll(e,alt.f,null.f,div=div)
lst <- test.results(res,e,pData(e)$treat,"U600","U200",div,
                    alpha=0.05,minSpC=2,minLFC=log2(1.8),
                    method="BH")

## ----volc1---------------------------------------------------------------
plot(lst$tres$LogFC, -log10(lst$tres$p.value))
plot(lst$tres$LogFC, -log10(lst$tres$p.value),
     xlim = c(-3, 3))
grid()

## ----volc2---------------------------------------------------------------
## Plot
res.volcanoplot(lst$tres,
                max.pval=0.05,
                min.LFC=1,
                maxx=3,
                maxy=NULL,
                ylbls=4)

## ----msmsedapca----------------------------------------------------------
library("msmsEDA")
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
lst <- counts.pca(msnset, wait=FALSE)

## ----pca-----------------------------------------------------------------
pcadata <- lst$pca$x[, 1:2]
head(pcadata)
plot(pcadata[, 1], pcadata[, 2],
     xlab = "PCA1", ylab = "PCA2")
grid()

## ----mkplottab, echo=FALSE-----------------------------------------------
plotfuns <- rbind(c("scatterplots", "plot", "xyplot", "geom_point"),
                  c("histograms", "hist", "histgram", "geom_histogram"),
                  c("density plots", "plot(density())", "densityplot", "geom_density"),
                  c("boxplots", "boxplot", "bwplot", "geom_boxplot"),
                  c("violin plots", "vioplot::vioplot", "bwplot(..., panel = panel.violin)", "geom_violin"),
                  c("line plots", "plot, matplot", "xyploy, parallelplot", "geom_line"),
                  c("bar plots", "barplot", "barchart", "geom_bar"),
                  c("pie charts", "pie", "", "geom_bar with polar coordinates"),
                  c("dot plots", "dotchart", "dotplot", "geom_point"),
                  c("stip plots", "stripchart", "stripplot", "goem_point"),
                  c("dendrogramms", "plot(hclust())", "latticeExtra package", "ggdendro package"),
                  c("heatmaps", "image, heatmap", "levelplot", "geom_tile"))


colnames(plotfuns) <- c("plot type", "traditional", "lattice", "ggplot2")

## ----, plottab-----------------------------------------------------------
kable(plotfuns)

## ----tandata-------------------------------------------------------------
library("pRolocdata")
data(tan2009r1)

## ----histex--------------------------------------------------------------
x <- exprs(tan2009r1)[, 1]

## ----histplot, fig.height=7, fig.width=14--------------------------------
par(mfrow = c(1, 2))
hist(x)
plot(density(x))

## ----bxplot, fig.width=7, fig.height=10----------------------------------
library("beanplot")
x <- exprs(tan2009r1)
par(mfrow = c(2, 1))
boxplot(x)
beanplot(x[, 1], x[, 2], x[, 3], x[, 4], log = "")

## ----matplotex, fig.width=7, fig.height=10-------------------------------
er <- fData(tan2009r1)$markers == "ER"
mt <- fData(tan2009r1)$markers == "mitochondrion"

par(mfrow = c(2, 1))
matplot(t(x[er, ]), type = "b", col = "red", pch = 1, lty = 1)
matplot(t(x[mt, ]), type = "b", col = "steelblue", pch = 1, lty = 1)

## ----mrktab--------------------------------------------------------------
x <- table(fData(tan2009r1)$markers)
x

## ----mrkplot, fig.height=7, fig.width=12---------------------------------
par(mfrow = c(1, 2))
barplot(x)
dotchart(x)

## ----hmap----------------------------------------------------------------
heatmap(exprs(tan2009r1))

## ----image, fig.width=14, fig.height=7-----------------------------------
par(mfrow = c(1, 2))
x <- matrix(1:9, ncol = 3)
image(x)
image(tan2009r1)

## ----dendro--------------------------------------------------------------
d <- dist(t(exprs(tan2009r1))) ## distance between samples
hc <- hclust(d) ## hierarchical clustering
plot(hc) ## visualisation

## ----mapsprep------------------------------------------------------------
library("mzR")
mzf <- pxget(px1, 6)
ms <- openMSfile(mzf)

hd <- header(ms)
ms1 <- which(hd$msLevel == 1)

rtsel <- hd$retentionTime[ms1] / 60 > 30 & hd$retentionTime[ms1] / 60 < 35
library("MSnbase")
(M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd))

## ----mapsheat------------------------------------------------------------
library("lattice")
ff <- colorRampPalette(c("yellow", "steelblue"))
trellis.par.set(regions=list(col=ff(100)))
plot(M, aspect = 1, allTicks = FALSE)

## ----maps3D--------------------------------------------------------------
M@map[msMap(M) == 0] <- NA
plot3D(M, rgl = FALSE)

## ----rglmap, eval=FALSE--------------------------------------------------
## library("rgl")
## plot3D(M, rgl = TRUE)

## ----msdetails-----------------------------------------------------------

lout <- matrix(NA, ncol = 10, nrow = 8)
lout[1:2, ] <- 1
for (ii in 3:4)
    lout[ii, ] <- c(2, 2, 2, 2, 2, 2, 3, 3, 3, 3)
lout[5, ] <- rep(4:8, each = 2)
lout[6, ] <- rep(4:8, each = 2)
lout[7, ] <- rep(9:13, each = 2)
lout[8, ] <- rep(9:13, each = 2)

i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
ms2 <- (i+1):(j-1)

layout(lout)

par(mar=c(4,2,1,1))
chromatogram(ms)
abline(v = hd[i, "retentionTime"], col = "red")


par(mar = c(3, 2, 1, 0))
plot(peaks(ms, i), type = "l", xlim = c(400, 1000))
legend("topright", bty = "n",
       legend = paste0(
           "Acquisition ", hd[i, "acquisitionNum"],  "\n",
           "Retention time ", formatRt(hd[i, "retentionTime"])))
abline(h = 0)
abline(v = hd[ms2, "precursorMZ"],
       col = c("#FF000080",
           rep("#12121280", 9)))

par(mar = c(3, 0.5, 1, 1))
plot(peaks(ms, i), type = "l", xlim = c(521, 522.5),
     yaxt = "n")
abline(h = 0)
abline(v = hd[ms2, "precursorMZ"], col = "#FF000080")

##par(mar = omar)
par(mar = c(2, 2, 0, 1))
for (ii in ms2) {
    p <- peaks(ms, ii)
    plot(p, xlab = "", ylab = "", type = "h", cex.axis = .6)
    legend("topright", legend = paste0("Prec M/Z\n",
                           round(hd[ii, "precursorMZ"], 2)),
           bty = "n", cex = .8)
}


## ----maps3D2-------------------------------------------------------------
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)

## ----barcode, fig.height=2, fig.width=12---------------------------------

par(mar=c(4,1,1,1))
image(t(matrix(hd$msLevel, 1, nrow(hd))),
      xlab="Retention time",
      xaxt="n", yaxt="n", col=c("black","steelblue"))
k <- round(range(hd$retentionTime) / 60)
nk <- 5
axis(side=1, at=seq(0,1,1/nk), labels=seq(k[1],k[2],k[2]/nk))


## ----anim1, eval=FALSE---------------------------------------------------
## library("animation")
## an1 <- function() {
##     for (i in seq(0, 5, 0.2)) {
##         rtsel <- hd$retentionTime[ms1] / 60 > (30 + i) &
##             hd$retentionTime[ms1] / 60 < (35 + i)
##         M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)
##         M@map[msMap(M) == 0] <- NA
##         print(plot3D(M, rgl = FALSE))
##     }
## }
## 
## saveGIF(an1(), movie.name = "msanim1.gif")

## ----anim2, eval=FALSE---------------------------------------------------
## an2 <- function() {
##     for (i in seq(0, 2.5, 0.1)) {
##         rtsel <- hd$retentionTime[ms1] / 60 > 30 & hd$retentionTime[ms1] / 60 < 35
##         mz1 <- 520 + i
##         mz2 <- 522 + i
##         M <- MSmap(ms, ms1[rtsel], mz1, mz2, .005, hd)
##         M@map[msMap(M) == 0] <- NA
##         print(plot3D(M, rgl = FALSE))
##     }
## }
## 
## saveGIF(an2(), movie.name = "msanim2.gif")

## ----msnbviz-------------------------------------------------------------
library("MSnbase")
data(itraqdata)
itraqdata2 <- pickPeaks(itraqdata, verbose = FALSE)
plot(itraqdata[[25]], full=TRUE, reporters = iTRAQ4)
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 1, 1))
plot(itraqdata2[[25]], itraqdata2[[28]], sequences = rep("IMIDLDGTENK", 2))

## ----protviz-------------------------------------------------------------
library("protViz")
data(msms)

fi <- fragmentIon("TAFDEAIAELDTLNEESYK")
fi.cyz <- as.data.frame(cbind(c=fi[[1]]$c, y=fi[[1]]$y, z=fi[[1]]$z))
     
p <- peakplot("TAFDEAIAELDTLNEESYK",
              spec = msms[[1]],
              fi = fi.cyz,
              itol = 0.6,
              ion.axes = FALSE)

## ----strp----------------------------------------------------------------
str(p)

## ----mqraw---------------------------------------------------------------
library("MALDIquant")

data("fiedler2009subset", package="MALDIquant")

plot(fiedler2009subset[[14]])

## ----mqestimatebaseline--------------------------------------------------
transformedSpectra <- transformIntensity(fiedler2009subset, method = "sqrt")
smoothedSpectra <- smoothIntensity(transformedSpectra, method = "SavitzkyGolay")

plot(smoothedSpectra[[14]])
lines(estimateBaseline(smoothedSpectra[[14]]), lwd = 2, col = "red")

## ----mqremovebaseline----------------------------------------------------
rbSpectra <- removeBaseline(smoothedSpectra)
plot(rbSpectra[[14]])

## ----mqpeaks-------------------------------------------------------------
cbSpectra <- calibrateIntensity(rbSpectra, method = "TIC")
peaks <- detectPeaks(cbSpectra, SNR = 5)

plot(cbSpectra[[14]])
points(peaks[[14]], col = "red", pch = 4, lwd = 2)

## ----mqlabelpeaks, echo = -(1:2)-----------------------------------------
plot(cbSpectra[[14]])
points(peaks[[14]], col = "red", pch = 4, lwd = 2)
top5 <- intensity(peaks[[14]]) %in% sort(intensity(peaks[[14]]),
                                         decreasing = TRUE)[1:5]
labelPeaks(peaks[[14]], index = top5, avoidOverlap = TRUE)

## ----mqwarp, fig.keep = "last"-------------------------------------------
par(mfrow = c(2, 2))
warpingFunctions <-
    determineWarpingFunctions(peaks,
                              tolerance = 0.001,
                              plot = TRUE,
                              plotInteractive = TRUE)
par(mfrow = c(1, 1))
warpedSpectra <- warpMassSpectra(cbSpectra, warpingFunctions)
warpedPeaks <- warpMassPeaks(peaks, warpingFunctions)

## ----mqwarped------------------------------------------------------------
sel <- c(2, 10, 14, 16)
xlim <- c(4180, 4240)
ylim <- c(0, 1.9e-3)
lty <- c(1, 4, 2, 6)

par(mfrow = c(1, 2))
plot(cbSpectra[[1]], xlim = xlim, ylim = ylim, type = "n")

for (i in seq(along = sel)) {
  lines(peaks[[sel[i]]], lty = lty[i], col = i)
  lines(cbSpectra[[sel[i]]], lty = lty[i], col = i)
}

plot(cbSpectra[[1]], xlim = xlim, ylim = ylim, type = "n")

for (i in seq(along = sel)) {
  lines(warpedPeaks[[sel[i]]], lty = lty[i], col = i)
  lines(warpedSpectra[[sel[i]]], lty = lty[i], col = i)
}
par(mfrow = c(1, 1))

## ----mqims, cache=FALSE, eval=FALSE, warning=FALSE-----------------------
## library("MALDIquant")
## library("MALDIquantForeign")
## 
## spectra <- importBrukerFlex("http://files.figshare.com/1106682/MouseKidney_IMS_testdata.zip", verbose = FALSE)
## 
## spectra <- smoothIntensity(spectra, "SavitzkyGolay",  halfWindowSize = 8)
## spectra <- removeBaseline(spectra, method = "TopHat", halfWindowSize = 16)
## spectra <- calibrateIntensity(spectra, method = "TIC")
## avgSpectrum <- averageMassSpectra(spectra)
## avgPeaks <- detectPeaks(avgSpectrum, SNR = 5)
## 
## avgPeaks <- avgPeaks[intensity(avgPeaks) > 0.0015]
## 
## oldPar <- par(no.readonly = TRUE)
## layout(matrix(c(1,1,1,2,3,4), nrow = 2, byrow = TRUE))
## plot(avgSpectrum, main = "mean spectrum",
##      xlim = c(3000, 6000), ylim = c(0, 0.007))
## lines(avgPeaks, col = "red")
## labelPeaks(avgPeaks, cex = 1)
## 
## par(mar = c(0.5, 0.5, 1.5, 0.5))
## for (i in seq(along = avgPeaks)) {
##   range <- mass(avgPeaks)[i] + c(-1, 1)
##   plotImsSlice(spectra, range = range,
##                main = paste(round(range, 2), collapse = " - "))
## }
## par(oldPar)

## ----ims-shiny, eval=FALSE-----------------------------------------------
## library("shiny")
## runGitHub("sgibb/ims-shiny")

## ----spatprot------------------------------------------------------------
library("pRoloc")
library("pRolocdata")

data(tan2009r1)

## these params use class weights
fn <- dir(system.file("extdata", package = "pRoloc"),
          full.names = TRUE, pattern = "params2.rda")
load(fn)

setStockcol(NULL)
setStockcol(paste0(getStockcol(), 90))

w <- table(fData(tan2009r1)[, "pd.markers"])
(w <- 1/w[names(w) != "unknown"])
tan2009r1 <- svmClassification(tan2009r1, params2,
                               class.weights = w,
                               fcol = "pd.markers")
ptsze <- exp(fData(tan2009r1)$svm.scores) - 1

## ----spatplot, fig.width=12, fig.height=6--------------------------------

lout <- matrix(c(1:4, rep(5, 4)), ncol = 4, nrow = 2)
layout(lout)
cls <- getStockcol()
par(mar = c(4, 4, 1, 1))
plotDist(tan2009r1[which(fData(tan2009r1)$PLSDA == "mitochondrion"), ],
         markers = featureNames(tan2009r1)[which(fData(tan2009r1)$markers.orig == "mitochondrion")],
         mcol = cls[5])
legend("topright", legend = "mitochondrion", bty = "n")
plotDist(tan2009r1[which(fData(tan2009r1)$PLSDA == "ER/Golgi"), ],
         markers = featureNames(tan2009r1)[which(fData(tan2009r1)$markers.orig == "ER")],
         mcol = cls[2])
legend("topright", legend = "ER", bty = "n")
plotDist(tan2009r1[which(fData(tan2009r1)$PLSDA == "ER/Golgi"), ],
         markers = featureNames(tan2009r1)[which(fData(tan2009r1)$markers.orig == "Golgi")],
         mcol = cls[3])
legend("topright", legend = "Golgi", bty = "n")
plotDist(tan2009r1[which(fData(tan2009r1)$PLSDA == "PM"), ],
         markers = featureNames(tan2009r1)[which(fData(tan2009r1)$markers.orig == "PM")],
         mcol = cls[8])
legend("topright", legend = "PM", bty = "n")
plot2D(tan2009r1, fcol = "svm", cex = ptsze, method = "kpca")
addLegend(tan2009r1, where = "bottomleft", fcol = "svm", bty = "n")


## ----si------------------------------------------------------------------
print(sessionInfo(), locale = FALSE)

