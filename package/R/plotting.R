#' Plots a density and true/predicted thresholds.
#'
#' @param densX Density x values.
#' @param densY Density y values.
#' @param thresh True threshold.
#' @param predicted Predicted threshold.
#'
#' @export
plotDensThresh <- function(densX, densY, thresh = NULL, predicted = NULL)
{
    plot(
         densX,
         densY,
         xlab = 'Channel',
         ylab = 'Density',
         type = 'o'
        );

    if(!is.null(thresh))
    {
        abline(v = thresh, col = 'red', lwd = 2, lty = 2)
    }    

    if(!is.null(predicted))
    {
        abline(v = predicted, col = 'blue', lwd = 2, lty = 2)
    }
}

#' Plots a gates as given by a TrainingGate object.
#'
#' @param trainGate The TrainingGate object to plot.
#' @param numFeatures Number of features to use.
#'
#' @export
plotGate <- function(trainGate, numFeatures = 512)
{
    par(mfrow = c(2, 1))

    neg <- trainGate@negate

    minX <- min(trainGate@densitiesA[[as.character(numFeatures)]]$x)
    maxX <- max(trainGate@densitiesA[[as.character(numFeatures)]]$x)
    minY <- min(trainGate@densitiesA[[as.character(numFeatures)]]$y)
    maxY <- max(trainGate@densitiesA[[as.character(numFeatures)]]$y)
    low <- if (!is.nan(trainGate@thresholdALow)) trainGate@thresholdALow else (if (!neg) minX else maxX)
    high <- if (!is.nan(trainGate@thresholdAHigh)) trainGate@thresholdAHigh else (if (!neg) maxX else minX)
    plot(
         trainGate@densitiesA[[as.character(numFeatures)]]$x,
         trainGate@densitiesA[[as.character(numFeatures)]]$y,
         xlab = 'Channel A',
         ylab = 'Density',
         type = 'o'
        );
    abline(v = low, col = 'red', lwd = 2, lty = 2)
    abline(v = high, col = 'red', lwd = 2, lty = 2)

    minX <- min(trainGate@densitiesB[[as.character(numFeatures)]]$x)
    maxX <- max(trainGate@densitiesB[[as.character(numFeatures)]]$x)
    minY <- min(trainGate@densitiesB[[as.character(numFeatures)]]$y)
    maxY <- max(trainGate@densitiesB[[as.character(numFeatures)]]$y)
    low <- if (!is.nan(trainGate@thresholdBLow)) trainGate@thresholdBLow else (if (!neg) minX else maxX)
    high <- if (!is.nan(trainGate@thresholdBHigh)) trainGate@thresholdBHigh else (if (!neg) maxX else minX)
    plot(
         trainGate@densitiesB[[as.character(numFeatures)]]$x,
         trainGate@densitiesB[[as.character(numFeatures)]]$y,
         xlab = 'Channel B',
         ylab = 'Density',
         type = 'o'
        );
    abline(v = low, col = 'red', lwd = 2, lty = 2)
    abline(v = high, col = 'red', lwd = 2, lty = 2)
}