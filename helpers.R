library(flowCore)

# A helper function for the C-style printf() function.
printf <- function(...) invisible(cat(sprintf(...)))

normalizePopulationName <- function(population)
# Normalizes a population name.
# First, all '+' and '-' signs are converted to 'p' and 'n', respectively.
# Second, All non-characters and non-numbers are removed.
# Third, everything is converted to lower case.
#
# Args:
#   population: The unnormalized population name string.
#
# Returns:
#   The normalized population string.
{
    populationNormalized <- str_replace_all(population, '\\+', 'p')
    populationNormalized <- str_replace_all(populationNormalized, '-', 'n')
    populationNormalized <- tolower(str_replace_all(populationNormalized, '[^a-zA-Z0-9]+', ''))
    populationNormalized
}

estimateDensity <- function(data, n)
# Estimates a smoothed density from a given vector of numbers.
#
# Args:
#   data: The numbers to calculated the density for.
#   n: Number of features for the density.
#
# Returns:
#   The estimated density.
{
  dens <- density(data[which(!is.na(data))], n = n)
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

rotateData <- function(data, chans=NULL, theta=NULL)
# Rotates FCS data by theta.
#
# Args:
#   data: The flowFrame to rotate or the ecpression matrix.
#   chans: Only rotate the given channels.
#   theta: Amount of rotation to apply.
#
# Returns:
#   A rotated version of the input.
{
    if (class(data)== "flowFrame" & !is.null(chans))
    {
        dataNew <- exprs(data)[,chans]
        if (is.null(theta))
        {
            regSlope <- atan(lm(dataNew[,1] ~ dataNew[,2])$coefficients[2])
            theta <- pi/2 - regSlope
        }
        dataNew <- dataNew %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
        exprs(data)[,chans] <- dataNew
    }else{
        data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    }
    return(list(data=data,theta=theta))
}

plotDensThresh <- function(densX, densY, thresh = NULL, predicted = NULL)
# Plots a density and true/predicted thresholds.
#
# Args:
#   densX: Density x values.
#   densY: Density y values.
#   thresh: True threshold.
#   predicted: Predicted threshold.
#
# Returns:
#   Nothing.
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

plotGate <- function(trainGate, numFeatures = 512)
# Plots a gates as given by a TrainingGate object.
#
# Args:
#   trainGate: The TrainingGate object to plot.
#   numFeatures: Number of features to use.
#
# Returns:
#   Nothing.
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

getMeanProportion <- function(tr, population)
# Reads proportions for a given population and return the mean proportion.
#
# Args:
#   tr: The LearningSet object to calculate the statistic for.
#   population: The population name.
#
# Returns:
#   The mean proportion.
{
    populationNormalized <- normalizePopulationName(population)

    n <- length(tr[[populationNormalized]]@samples)

    cl <- makeCluster(detectCores(), type = "FORK")

    tryCatch({
        props <- t(parSapply(cl, 1:n, function(i) 
        {
            tg <- readRDS(paste0('trainingFiles/', tr[[populationNormalized]]@samples[[i]]))

            gateAssignments <- tg[[population]]@gateAssignments     

            sum(gateAssignments) / length(gateAssignments)
        }))

        return(mean(props))
    }, error = function(e) {
        print(e)
    }, finally = {
        stopCluster(cl) 
    }) 
}