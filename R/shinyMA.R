##' @import shiny
shinyMA <- function() {
    qnt <- NULL
    data(qnt, envir = environment())
    findClosest <- function(input) {
        if (is.null(input$plotma_click)) {
            1
        } else {
            pt <- c(input$plotma_click$x, y=input$plotma_click$y)
            sqdists <- colMeans((t(fData(qnt)[, c("baseMean", "log2FC")]) - pt)^2)
            which.min(sqdists)
        }
    }

    app <- list(
        ui = fluidPage(
            titlePanel("MA-plot + expression plot with shiny"),

            ## flow layout fills out left to right then down
            flowLayout(
                plotOutput("plotma", clickId="plotma_click", width=400, height=400),
                plotOutput("plotcounts", width=400, height=400)
            )
        ),

        server = function(input, output) {
            ## MA-plot
            output$plotma <- renderPlot({
                ## process possible input
                idx <- findClosest( input )
                par(mar = c(5,5,3,2), cex.lab = 1.5)
                ## MA-plot of all genes
                plot(fData(qnt)$baseMean,
                     fData(qnt)$log2FC,
                     xlab = "A", ylab = "M")
                grid()
                abline(v = 0)
                ## add circle for the selected point
                points(fData(qnt)$baseMean[idx],
                       fData(qnt)$log2FC[idx],
                       col="dodgerblue", cex=3, lwd=3)
            })

            ## counts plot
            output$plotcounts <- renderPlot({
                ## process possible input
                idx <- findClosest( input )
                par(mar = c(5,5,3,2), cex.lab = 1.5)
                ## plot the counts for the selected gene
                i <- exprs(qnt)[idx, ]
                plot(1:6, i,
                     col = c("red", "steelblue")[as.numeric(qnt$conditions)],
                     pch = 19,
                     main = paste0(fData(qnt)$accession[idx],
                         " (", fData(qnt)$sequence[idx], ")"),
                     xlab = "samples",
                     ylab = "intensity")
                abline(h = mean(i[qnt$conditions == "A"]),
                       col = "red", lty = "dotted")
                abline(h = mean(i[qnt$conditions == "B"]),
                       col = "steelblue", lty = "dotted")
            })
        }
    )
    runApp(app)
}
