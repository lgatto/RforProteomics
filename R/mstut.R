##' This is the constructor function to generate a set of ions that
##' can later be analysed with `analyse()` and detected with
##' `detect()`. 
##'
##' @title Create, analyse and detect ions
##' @param npeaks A `numeric` scalar defining the number of unique
##'     peaks (M/Z values). Default is 10.
##' @param mzrange A `numeric` of length 2 defining the range of
##'     possible M/Z values. Default is `c(100, 1000)`.
##' @param nimg A `numeric` scalar. When analysing the ions, their
##'     separation along their M/Z values will be split along a
##'     sequence of length `nimg`. Default is 100.
##' @return An object of class `ions`.
##' @author Laurent Gatto
##' @rdname mstut
##' @aliases mstut
##' @export
##' @examples
##' set.seed(1L)
##' x <- new_ions(nimg = 5)
##' x
##' analyse(x)
##' detect(x)
##' spectrum(x)
new_ions <- function(npeaks = 10,
                     mzrange = c(100, 1000),
                     nimg = 100) {
    ## peaks
    mzs <- stats::runif(npeaks, min = min(mzrange), max = max(mzrange))    
    maxint <- 10
    k <- sample(maxint, npeaks, replace = TRUE)    
    mzs <- sample(rep(mzs, k))
    N <- length(mzs)

    ## visuals
    cex <- mzs / max(mzs) * 2

    ## ms data
    ys <- seq(0.05, 0.95, length = N)
    x0 <- rep(100, N) ## start
    msdata <- mapply(seq, x0, mzs, length = nimg)

    ## spectrum
    smzs <- sort(mzs)
    sp <- data.frame(MZ = unique(smzs),
                     Intensity = as.vector(table(smzs))/maxint)

    structure(list(mzrange = range(mzs),
                   ys = ys,
                   msdata = msdata,
                   size = cex,
                   N = N,
                   spectrum = sp),
              class = "ions")
}


##' @export
print.ions <- function(x, ...) {
    cat("Object of class 'ions':\n")
    cat("# analyze(.); detect(.); spectrum(.)\n")
}

##' @param x An object of class `ions`.
##' @param sleep How much time to wait before producing the next plot.
##' @return `analyse`, `detect` and `spectrum` are used for their side
##'     effect or producing plots. They all invisibly return `NULL`.
##' @export
##' @importFrom graphics plot
##' @rdname mstut
analyse <- function(x, sleep = 0.1) {
    stopifnot(inherits(x, "ions"))
    n <- nrow(x$msdata)
    sapply(seq_len(nrow(x$msdata)),
          function(i) {
              plot(x$msdata[i, ], x$ys, 
                   xlim = x$mzrange, ylim = c(0, 1),
                   yaxt = "n", xlab = "M/Z", ylab = "Analytes",
                   main = paste0("Analyser (", i, "/", n, ")"),
                   cex = x$size)
              Sys.sleep(sleep)
          })
    invisible(NULL)
}

##' @rdname mstut
##' @export
analyze <- analyse

##' @param new A `logical` scalar, indicating if the separated ions
##'     (last frame of calling `analyse) should be plotting, or
##'     whether the detection should be overlaid. Default is `FALSE`,
##'     to add the plot on top of the opened device.
##' @export
##' @importFrom graphics grid lines
##' @rdname mstut
detect <- function(x, new = FALSE) {
    stopifnot(inherits(x, "ions"))
    if (new)
        plot(x$msdata[nrow(x$msdata),], x$ys, 
             xlim = x$mzrange, ylim = c(0, 1),
             yaxt = "n", xlab = "M/Z", ylab = "Analytes",
             main = "Analyser", cex = x$size)        
    grid()
    lines(x$spectrum$MZ, x$spectrum$Intensity,
          col = "red", type = "h", lwd = 2) 

}

##' @param ... Additional arguments passed to [graphics::plot()].
##' @export
##' @rdname mstut
spectrum <- function(x, ...) {
    plot(x$spectrum, type = "h",
         ylim = c(0, 1),
         lwd = 2, ...)
    grid()
}

## empty_ms <- function(main = "Analyser") {
##     plot(NA, type = "n", xlim = c(100, 1000), ylim = c(0, 1),
##          yaxt = "n", xlab = "M/Z", ylab = "Analytes",
##          main = main)
## }





