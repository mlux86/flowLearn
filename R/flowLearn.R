setClass("GatingInfo", representation(population = "character", parent = "character", channels = "numeric")) # defined for each panel

#' An S4 class to represent density data used by flowLearn.
#'
#' @slot data A data frame where each row represents a density of a particular channel, together with additional information.
#'
#'       Columns of the data frame are:
#'
#'       - fcs: FCS file of the density.
#'
#'       - population: Analyzed population of the density.
#'
#'       - channelIdx: Identifier of the channel for which the density was calculated.
#'
#'       - <numFeatures> x-values of the density
#'
#'       - <numFeatures> y-values of the density
#'
#'       - gate.low: Lower threshold for this channel (defaults to NaN)
#'
#'       - gate.high: Upper threshold for this channel (defaults to NaN)
#'
#' @slot numFeatures Number of used features.
#'
#' @examples
#' str(flSampleDensdat) # summarizes the sample DensityData object provided by flowLearn
#'
#' @export
DensityData <- setClass("DensityData", slots = c(data = "data.frame", numFeatures = 'numeric'), prototype = list(data = data.frame(), numFeatures = 512))

#' Getter for the "data" slot
#'
#' @param obj A DensityData object.
#'
#' @return The internal data frame "data"
#'
#' @examples
#' flData(flSampleDensdat)
#'
#' @export
setGeneric(name = "flData", def = function(obj) { standardGeneric("flData") })
setMethod(f = "flData", signature = "DensityData", definition = function(obj) { return(obj@data) })

#' Getter for the "numFeatures" slot
#'
#' @param obj A DensityData object.
#'
#' @return The number of features used by this DensityData
#'
#' @examples
#' flNumFeatures(flSampleDensdat)
#'
#' @export
setGeneric(name = "flNumFeatures", def = function(obj) { standardGeneric("flNumFeatures") })
setMethod(f = "flNumFeatures", signature = "DensityData", definition = function(obj) { return(obj@numFeatures) })

#' Initializes a DensityData object with default values
#'
#' @param obj An uninitialized DensityData object.
#'
#' @return The initialized DensityData Object.
#'
#' @examples
#' str(flInit(new('DensityData'))) # summarizes an empty DensityData object that has just been initialized
#'
#' @export
setGeneric(name = "flInit", def = function(obj) { standardGeneric("flInit") })

setMethod(f = "flInit", signature = "DensityData",
                      definition = function(obj)
                      {
                   	  	  df1 <- data.frame(fcs = character(0), population = character(0), channelIdx = integer(0), stringsAsFactors = FALSE)
                   	  	  df2 <- setNames(replicate(obj@numFeatures, integer(0), simplify = FALSE), sapply(1:obj@numFeatures, function(i) { paste0('x.', i) }))
                   	  	  df3 <- setNames(replicate(obj@numFeatures, integer(0), simplify = FALSE), sapply(1:obj@numFeatures, function(i) { paste0('y.', i) }))
                   	  	  df4 <- data.frame(gate.low = integer(0), gate.high = integer(0), stringsAsFactors = FALSE)
                   	  	  obj@data <- cbind(df1, df2, df3, df4)
                          return(obj)
                      }
                      )

#' Adds one row (i.e. one density) to an existing DensityData object.
#'
#' @param obj The existing DensityData object.
#' @param fcs The fcs file name of the new density.
#' @param population The population of the new density.
#' @param channelIdx The channel index of the new density.
#' @param parentDensX The parent population's density x-values of the gated population.
#' @param parentDensY The parent population's density y-values of the gated population.
#' @param gate.low The lower threshold of the gated channel (optional, default: NaN).
#' @param gate.high The upper threshold of the gated channel (optional, default: NaN).
#'
#' @return The supplied DensityData obj, including a new row with the supplied parameters.
#'
#' @examples
#' densdat <- flInit(new('DensityData')) # create empty DensityData object
#' densdat <- flAdd(densdat, 'mine.fcs', 'CD43+', 1, 1:512, runif(n = 512, min = 1, max = 10), 2, 7) # add one channel density
#'
#' @export
setGeneric(name = "flAdd", def = function(obj, fcs, population, channelIdx, parentDensX, parentDensY, gate.low = NaN, gate.high = NaN) { standardGeneric("flAdd") })

setMethod(f = "flAdd", signature = "DensityData",
                     definition = function(obj, fcs, population, channelIdx, parentDensX, parentDensY, gate.low = NaN, gate.high = NaN)
                     {
                   	  	  df1 <- data.frame(fcs = fcs, population = population, channelIdx = channelIdx, stringsAsFactors = FALSE)
                   	  	  df2 <- as.data.frame(matrix(parentDensX, nrow = 1), stringsAsFactors = FALSE); colnames(df2) <- sapply(1:obj@numFeatures, function(i) { paste0('x.', i) })
                   	  	  df3 <- as.data.frame(matrix(parentDensY, nrow = 1), stringsAsFactors = FALSE); colnames(df3) <- sapply(1:obj@numFeatures, function(i) { paste0('y.', i) })
                   	  	  df4 <- data.frame(gate.low = gate.low, gate.high = gate.high, stringsAsFactors = FALSE)
                   	  	  obj@data <- rbind(obj@data, cbind(df1, df2, df3, df4))
                          return(obj)
                     }
                     )

#' Adds densities from a flowFrame that represents a parent cell population.
#'
#' @param obj A DensityData object.
#' @param parentFlowFrame An object of type flowFrame representing the parent population.
#' @param channelIndices A vector of indices of channels to use for gating the child population.
#' @param populationName The name of the child population.
#' @param fcsName The different name for the FCS file. By default parentFlowFrame@description$"$FIL" is used.
#'
#' @return The new DensityData object containing one density for each channel specified.
#'
#' @examples
#' densdat <- flInit(new('DensityData'))
#' flSampleFlowFrame <- readRDS(gzcon(url('https://raw.githubusercontent.com/mlux86/flowLearn/master/extra/data/flSampleFlowFrame.rds')))
#' densdat <- flAddFlowFrame(densdat, flSampleFlowFrame, c(1,2), 'Granulocyte Pre')
#'
#' @export
setGeneric(name = "flAddFlowFrame", def = function(obj, parentFlowFrame, channelIndices, populationName, fcsName = NULL) { standardGeneric("flAddFlowFrame") })

setMethod(f = "flAddFlowFrame", signature = "DensityData",
                      definition = function(obj, parentFlowFrame, channelIndices, populationName, fcsName = NULL)
                      {
                          if (is.null(fcsName))
                          {
                              fcsName <- parentFlowFrame@description$"$FIL"
                          }

                          for (c in channelIndices)
                          {
                              dens <- flEstimateDensity(parentFlowFrame@exprs[, c], obj@numFeatures)
                              obj <- flAdd(obj, fcsName, populationName, c, dens$x, dens$y)
                          }

                          return(obj)
                      }
                      )

#' Returns a list(x,y) with density x and y values for the given DensityData object. Usually this is called on DensityData objects with only one row, to extract a density for one specific channel.
#'
#' @param obj The DensityData object.
#'
#' @return A list(x = matrix, y = matrix) with matrices x and y that represent x and y parts of densities for each row in obj.
#'
#' @examples
#' dens <- flGetDensity(  flAt(flSampleDensdat, 42)  ) # grab density at row index 42
#' flPlotDensThresh(dens)
#'
#' @export
setGeneric(name = "flGetDensity", def = function(obj) { standardGeneric("flGetDensity") })

setMethod(f = "flGetDensity", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  x <- as.matrix(obj@data[,4:(3+obj@numFeatures)])
                          	  y <- as.matrix(obj@data[,516:(515+obj@numFeatures)])
                          	  list(x = x, y = y)
                          }
                          )

#' Returns a matrix with two columns which represent lower and upper thresholds of the given DensityData object.
#'
#' @param obj The DensityData object.
#'
#' @return A matrix with columns "gate.low" and "gate.high", representing lower and upper thresholds on all rows/channels present obj.
#'
#' @examples
#' dd <- flAt(flSampleDensdat, 42) # grab density at row index 42
#' dens <- flGetDensity(dd) #
#' gt <- flGetGate(dd) # returns vector with "gate.low" and "gate.high"
#' flPlotDensThresh(dens, gt)
#'
#' @export
setGeneric(name = "flGetGate", def = function(obj) { standardGeneric("flGetGate") })

setMethod(f = "flGetGate", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  as.matrix(obj@data[,c('gate.low', 'gate.high')])
                          }
                          )

#' Returns a new DensityData object with a subset of entries of the supplied obj.
#'
#' @param fcs The FCS file to filter for.
#' @param population The population to filter for.
#' @param channel The channel to filter for.
#'
#' @return A new DensityData object with a subset of entries, defined by mysubset.
#'
#' @examples
#' # Get all densities for Population "hfa" and the first channel
#' dd <- flFind(flSampleDensdat, population = 'hfa', channelIdx = 1)
#'
#' @export
setGeneric(name = "flFind", def = function(obj, fcs, population, channelIdx) { standardGeneric("flFind") })

setMethod(f = "flFind", signature = "DensityData",
                          definition = function(obj, fcs, population, channelIdx)
                          {
                              idx <- replicate(nrow(obj@data), TRUE)
                              if (!missing(fcs))
                              {
                                  idx <- idx & obj@data$fcs == fcs
                              }
                              if (!missing(population))
                              {
                                  idx <- idx & obj@data$population == population
                              }
                              if (!missing(channelIdx))
                              {
                                  idx <- idx & obj@data$channelIdx == channelIdx
                              }
                              obj@data <- obj@data[idx,]
                          	  return(obj)
                          }
                          )

#' Returns a single density from the data slot at the specified row index.
#'
#' @param obj The DensityData object to be filtered.
#' @param idx The row index of the internal "data" data frame.
#'
#' @return A new DensityData object with only one row, given by the specified index.
#'
#' @examples
#' dd <- flAt(flSampleDensdat, 42) # grab 42nd entry of flSampleDensdat
#'
#' @export
setGeneric(name = "flAt", def = function(obj, idx) { standardGeneric("flAt") })

setMethod(f = "flAt", signature = "DensityData",
                          definition = function(obj, idx)
                          {
                              newData <- obj@data[idx, ]
                              rownames(newData) <- NULL
                          	  initialize(obj, data = newData)
                          }
                          )

#' Returns the size of a given DensityData object.
#'
#' @param obj The DensityData object.
#'
#' @return The size of obj, i.e. the number of rows of the internal "data" data frame.
#'
#' @examples
#' print(flSize(flSampleDensdat)) # print size of flSampleDensdat
#'
#' @export
setGeneric(name = "flSize", def = function(obj) { standardGeneric("flSize") })

setMethod(f = "flSize", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  nrow(obj@data)
                          }
                          )

#' Concatenates two DensityData objects.
#'
#' @param obj1 The first DensityData object.
#' @param obj2 The second DensityData object.
#'
#' @return A new DensityData object consisting of the internal "data" data frames of obj1 and obj2.
#'
#' @examples
#' print(flSize(flSampleDensdat))
#' dd <- flConcat(flSampleDensdat, flSampleDensdat)
#' print(flSize(dd))
#'
#' @export
setGeneric(name = "flConcat", def = function(obj1, obj2) { standardGeneric("flConcat") })

setMethod(f = "flConcat", signature = "DensityData",
                          definition = function(obj1, obj2)
                          {
                              DensityData(data = rbind(obj1@data, obj2@data), numFeatures = obj1@numFeatures)
                          }
                          )

#' Normalizes population names, such that they can be used, for example in file names.
#'
#' In particular, '+' is replaced by 'p' (positive), '-' is replaced by 'n' (negative), and all special characters are removed.
#'
#' @param population The (unnormalized) name of a population.
#'
#' @return The normalized population name.
#'
#' @export
#'
#' @examples
#' popName <- 'IGD+ CD27-'
#' popNameNormalized <- flNormalizePopulationName(popName)
#' print(popNameNormalized)
#'
flNormalizePopulationName <- function(population)
{
    populationNormalized <- gsub('\\+', 'p', population)
    populationNormalized <- gsub('-', 'n', populationNormalized)
    populationNormalized <- tolower(gsub('[^a-zA-Z0-9]+', '', populationNormalized))
    return(populationNormalized)
}

#' Estimates a density on a given data vector.
#'
#' This method calculates densities used by flowLearn. It uses the default stats::density function and smoothens the result using stats::spline.smooth.
#'
#' @param data The data vector for which the density is calculated.
#' @param n The number of density features.
#'
#' @return The calculated density object.
#'
#' @export
#'
#' @examples
#' x <- seq(-pi, pi, by = 0.01)
#' y <- sin(x)
#' noise <- runif(n = length(y), min = 0, max = 1)
#' noisy_y <- y + noise
#' par(mfrow = c(2,1))
#' plot(noisy_y)
#' plot(flEstimateDensity(noisy_y, 512))
#'
flEstimateDensity <- function(data, n)
{
  dens <- stats::density(data[which(!is.na(data))], n = n, from = 0)
  dens <- stats::smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

#' Given two R density objects, calculates a distance matrix based on the derivative DTW distance.
#'
#' Given two R density objects, this method uses the Derivative Dynamic Time Warping distance to calculate a distance matrix which can be used for alignment with DTW.
#'
#' @param densA The first density object, such as from R's density function.
#' @param densB The second density object, such as from R's density function.
#'
#' @return The Derivative Dynamic Time Warping distance matrix between densA and densB.
#'
#' @importFrom proxy dist
#' @export
#'
#' @examples
#' densA <- flGetDensity(flAt(flSampleDensdat, 42))
#' densB <- flGetDensity(flAt(flSampleDensdat, 43))
#' str(flDerivativeDtwDistanceMatrix(densA, densB))
#'
flDerivativeDtwDistanceMatrix <- function(densA, densB)
{
	lenA <- length(densA$y)
	lenB <- length(densB$y)

    if(lenA != lenB)
    {
        stop('Densities A and B must have same length.')
    }

    difA <- sapply(2:(lenA-1), function(i) (densA$y[i] - densA$y[i-1]) + ((densA$y[i+1] - densA$y[i-1]) / 2) / 2)
    difA <- c(difA[1], difA, difA[lenA-2])

    difB <- sapply(2:(lenB-1), function(i) (densB$y[i] - densB$y[i-1]) + ((densB$y[i+1] - densB$y[i-1]) / 2) / 2)
    difB <- c(difB[1], difB, difB[lenB-2])

    proxy::dist(difA, difB)
}

#' Wrapper function to carry out DDTW for flowLearn.
#'
#' Wrapper function to carry out Derivative Dynamic Time Warping for flowLearn. It uses Derivative DTW with a DTW specific step pattern, which is "typeIds" from the dtw package.
#'
#' @param densA The first density object, such as from R's density function.
#' @param densB The second density object, such as from R's density function.
#' @param ... Optional parameters given to dtw::dtw
#'
#' @return The Derivative Dynamic Time Warping object from the dtw package.
#'
#' @importFrom dtw dtw
#' @importFrom dtw typeIds
#' @export
#'
#' @examples
#' dd <- flFind(flSampleDensdat, population = 'cd3tcell', channelIdx = 1)
#' dtwObj <- flDtwMain(flGetDensity(flAt(dd, 1)), flGetDensity(flAt(dd, 2)))
#' plot(dtwObj)
#'
flDtwMain <- function(densA, densB, ...)
{
    ddtw <- dtw::dtw(flDerivativeDtwDistanceMatrix(densA, densB), step.pattern = dtw::typeIds, ...)

    ddtw$query <- as.double(densA$y)
    ddtw$reference <- as.double(densB$y)

    ddtw
}

#' Given a DensityData object, selects a number of prototype densities from it.
#'
#' Using PAM clustering on L1 density distances, this method returns a list of prototype indices which are rows in the given DensityData object.
#'
#' @param densdat The DensityData object. It should contain densities from the same channel only.
#' @param k The number of prototypes.
#'
#' @return Vector of indices in the densdat data that were determined to be prototypes.
#'
#' @importFrom proxy dist
#' @export
#'
#' @examples
#' dd <- flFind(flSampleDensdat, population = 'notplasma', channelIdx = 2)
#' protoIdx <- flSelectPrototypes(dd, 1)
#' print(protoIdx)
#'
flSelectPrototypes <- function(densdat, k)
{
    D <- proxy::dist(flGetDensity(densdat)$y, method = 'manhattan')

    protoIdx <- cluster::pam(D, k = k)$medoids

    return(protoIdx)
}

#' Aligns two densities A and B with each other and transfers a given reference threshold in B to A.
#'
#' Uses Derivative Dynamic Time Warping to align a given density A with a reference density B and transfers a reference threshold fromB to A.
#'
#' @param dens A R density object with unknown threshold.
#' @param refDens A reference R density object with known threshold.
#' @param refThreshold The threshold belonging to refDens.
#'
#' @return A double value representing the threshold that resulted from aligning both densities and transferring the reference threshold. If refThreshold is NA, this method returns NA as well.
#'
#' @export
#'
#' @examples
#' dd <- flFind(flSampleDensdat, population = 'plasma', channelIdx = 1)
#' dens <- flGetDensity(flAt(dd, 1))
#' trueThresh <- flGetGate(flAt(dd, 1))[1]
#' refDens <- flGetDensity(flAt(dd, 2))
#' refThresh <- flGetGate(flAt(dd, 2))[1]
#'
#' predictedThreshold <- flAlignThreshold(dens, refDens, refThresh)
#' flPlotDensThresh(dens, refThresh, trueThresh)
#'
flAlignThreshold <- function(dens, refDens, refThreshold)
{
    if (is.na(refThreshold))
    {
        return(NA)
    }

    dtwObj <- flDtwMain(dens, refDens)

    # find x index nearest to threshold in reference

    nearestRefXIdx <- which.min(abs( refDens$x - refThreshold ))

    # identify mapped indexes on aligned density

    idx <- dtwObj$index1[dtwObj$index2 == nearestRefXIdx]

    if (length(idx) > 1)
    {
        thresh <- mean(dens$x[idx])
    } else
    {
        thresh <- dens$x[idx]
    }

    thresh
}

#' Predicts thresholds for all densities in a DensityData object, using a set of prototypes.
#'
#' This method uses the prototypes specified in protoIdx to predict thresholds of all other (non-prototype indices) in densdat.
#' For each density, it determines the nearest prototype density and calls flAlignThreshold to transfer lower and upper thresholds.
#' Densities in densdat that are prototypes (specified by protoIdx) must have thresholds set.
#' Existing non-prototype thresholds in densdat are replaced by predicted thresholds.
#'
#' @param densdat A DensityData object containing prototypes and non-prototypes. Prototypes must have thresholds set. All densities in densdat have to be from the same channel.
#' @param protoIdx A set of indices that specify prototype rows in densdat.
#'
#' @return A copy of densdat where non-prototype densities have predicted thresholds set.
#'
#' @export
#'
#' @examples
#' dd <- flFind(flSampleDensdat, population = 'notplasma', channelIdx = 2)
#' protoIdx <- flSelectPrototypes(dd, 1)
#' ddp <- flPredictThresholds(dd, protoIdx)
#'
#' par(mfrow = c(3,1))
#' for(i in 1:3) flPlotDensThresh(flGetDensity(flAt(ddp, i)), flGetGate(flAt(dd, i)), flGetGate(flAt(ddp, i)))
#'
flPredictThresholds <- function(densdat, protoIdx)
{
	n <- flSize(densdat)
	nProto <- length(protoIdx)

	testIdx <- setdiff(1:n, protoIdx)
	nTest <- length(testIdx)

	D <- as.matrix(dist(flGetDensity(densdat)$y, method = 'manhattan'))

	cl <- parallel::makeCluster(parallel::detectCores(), type = "FORK")


  tryCatch({

      predicted <- t(parallel::parSapply(cl, testIdx, function(i)
      {
          j <- order(D[protoIdx, i])[1]
          g.l <- flAlignThreshold(flGetDensity(flAt(densdat, i)), flGetDensity(flAt(densdat, protoIdx[j])), flGetGate(flAt(densdat, protoIdx[j]))[1])
          g.h <- flAlignThreshold(flGetDensity(flAt(densdat, i)), flGetDensity(flAt(densdat, protoIdx[j])), flGetGate(flAt(densdat, protoIdx[j]))[2])
          c(g.l, g.h)
      }))

	    ndim <- ncol(densdat@data)
	    densdat@data[testIdx, c(ndim-1, ndim)] <- predicted
		return(densdat)

	}, error = function(e) {
		print(e)
	}, finally = {
		parallel::stopCluster(cl)
	})

}

#' Plots a given density and true/predicted thresholds, if supplied.
#'
#' Plots a given density and true/predicted thresholds, if supplied. True thresholds are drawn as red lines. Predicted thresholds are drawn as blue lines.
#' If no thresholds are specified, they are not drawn.
#'
#' @param dens A R density object such as returned from the stats::density function.
#' @param thresh The true threshold in the given density. (optional)
#' @param predicted The predicted threshold in the given density.
#' @param xlab The label of the x-axis (defaults to "Channel")
#' @param ylab The label of the y-axis (defaults to "Density")
#'
#' @return Nothing
#'
#' @export
#'
#' @examples
#' dd <- flFind(flSampleDensdat, population = 'myeloid', channelIdx = 1)
#' protoIdx <- flSelectPrototypes(dd, 1)
#' ddp <- flPredictThresholds(dd, protoIdx)
#'
#' par(mfrow = c(3,1))
#' for(i in 1:3) flPlotDensThresh(flGetDensity(flAt(ddp, i)), flGetGate(flAt(dd, i)), flGetGate(flAt(ddp, i)))
#'
flPlotDensThresh <- function(dens, thresh = NULL, predicted = NULL, xlab = "Channel", ylab = "Density")
{
    plot(
         dens$x,
         dens$y,
         xlab = xlab,
         ylab = ylab,
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

