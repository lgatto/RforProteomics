\name{qnt}
\alias{qnt}
\docType{data}

\title{
  PXD000001 example MSnSet
}

\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}

\usage{data("qnt")}

\format{
  An instance of class \code{MSnSet}
}

\details{

  In this TMT 6-plex experiment, four exogenous proteins were spiked
  into an equimolar \emph{Erwinia carotovora} lysate with varying
  proportions in each channel of quantitation; yeast enolase (ENO) at
  10:5:2.5:1:2.5:10, bovine serum albumin (BSA) at 1:2.5:5:10:5:1,
  rabbit glycogen phosphorylase (PHO) at 2:2:2:2:1:1 and bovin
  cytochrome C (CYT) at 1:1:1:1:1:2. Proteins were then digested,
  differentially labelled with TMT reagents, fractionated by reverse
  phase nanoflow UPLC (nanoACQUITY, Waters), and analysed on an LTQ
  Orbitrap Velos mass spectrometer (Thermo Scientific).  Files in
  multiple format will be used to illustrate the input/output
  capabilities that are available to the proteomics audience. The
  companion package provides dedicated functions to directly download
  the data.

  The data has been downloaded from the ProteomeXchange repository and
  imported into R as illustrated in the example.
}


\references{

  Laurent Gatto (2014). RforProteomcs: Companion package to the 'Using R
  and Bioconductor for proteomics data analysis' publicationR package
  version 1.3.1.

  Gatto L, Christoforou A. Using R and Bioconductor for proteomics data
  analysis. Biochim Biophys Acta. 2013 May 18. doi:pii:
  S1570-9639(13)00186-6. 10.1016/j.bbapap.2013.04.032. [Epub ahead of
  print] PubMed PMID: 23692960.

}

\examples{

\dontrun{
    library("rpx")
    px1 <- PXDataset("PXD000001")
    mztab <- pxget(px1, "PXD000001_mztab.txt")
    library("MSnbase")
    qnt <- readMzTabData(mztab, what = "PEP")
    sampleNames(qnt) <- reporterNames(TMT6)
    qnt$conditions <- factor(c("A", "A", "B", "B", "B", "A"))
    qnt <- filterNA(qnt)

    selA <- qnt$conditions == "A"

    fData(qnt)$log2FC <-
        log(rowMeans(exprs(qnt)[, selA]), 2) -
            log(rowMeans(exprs(qnt)[, !selA]), 2)
    fData(qnt)$baseMean <- log(rowMeans(exprs(qnt)), 10)
}

library("RforProteomics")
data(qnt)
head(exprs(qnt))
head(fData(qnt))
}

\keyword{datasets}
