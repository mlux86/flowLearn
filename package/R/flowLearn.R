library(dtw)
library(stringr)
library(parallel)
library(cluster)

setClass("GatingInfo", representation(population = "character", parent = "character", channels = "numeric")) # defined for each panel

#' @export
DensityData <- setClass("DensityData", slots = c(data = "data.frame", numFeatures = 'numeric'), prototype = list(data = data.frame(), numFeatures = 512))

setGeneric(name = "flInit", def = function(obj) { standardGeneric("flInit") })

setMethod(f = "flInit", signature = "DensityData",
                      definition = function(obj)
                      {
                   	  	  df1 <- data.frame(fcs = character(0), population = character(0), channelIdx = integer(0), stringsAsFactors = F)
                   	  	  df2 <- setNames(replicate(obj@numFeatures, integer(0), simplify = F), sapply(1:obj@numFeatures, function(i) { paste0('x.', i) }))
                   	  	  df3 <- setNames(replicate(obj@numFeatures, integer(0), simplify = F), sapply(1:obj@numFeatures, function(i) { paste0('y.', i) }))
                   	  	  df4 <- data.frame(gate.low = integer(0), gate.high = integer(0), stringsAsFactors = F)
                   	  	  obj@data <- cbind(df1, df2, df3, df4)
                          return(obj)
                      }
                      )

setGeneric(name = "flAdd", def = function(obj, fcs, population, channelIdx, parentDensX, parentDensY, gate.low = NaN, gate.high = NaN) { standardGeneric("flAdd") })

setMethod(f = "flAdd", signature = "DensityData",
                     definition = function(obj, fcs, population, channelIdx, parentDensX, parentDensY, gate.low = NaN, gate.high = NaN)
                     {
                   	  	  df1 <- data.frame(fcs = fcs, population = population, channelIdx = channelIdx, stringsAsFactors = F)
                   	  	  df2 <- as.data.frame(matrix(parentDensX, nrow = 1), stringsAsFactors = F); colnames(df2) <- sapply(1:obj@numFeatures, function(i) { paste0('x.', i) })
                   	  	  df3 <- as.data.frame(matrix(parentDensY, nrow = 1), stringsAsFactors = F); colnames(df3) <- sapply(1:obj@numFeatures, function(i) { paste0('y.', i) })
                   	  	  df4 <- data.frame(gate.low = gate.low, gate.high = gate.high, stringsAsFactors = F)
                   	  	  obj@data <- rbind(obj@data, cbind(df1, df2, df3, df4))
                          return(obj)
                     }
                     )

setGeneric(name = "flGetDensity", def = function(obj) { standardGeneric("flGetDensity") })

setMethod(f = "flGetDensity", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  x <- as.matrix(obj@data[,4:(3+obj@numFeatures)])
                          	  y <- as.matrix(obj@data[,516:(515+obj@numFeatures)])
                          	  list(x = x, y = y)
                          }
                          )

setGeneric(name = "flGetGate", def = function(obj) { standardGeneric("flGetGate") })

setMethod(f = "flGetGate", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  as.matrix(obj@data[,c('gate.low', 'gate.high')])
                          }
                          )

setGeneric(name = "flFind", def = function(obj, mysubset) { standardGeneric("flFind") })

setMethod(f = "flFind", signature = "DensityData",
                          definition = function(obj, mysubset)
                          {
                              obj@data <- subset(obj@data, eval(parse(text = mysubset)))
                              rownames(obj@data) <- NULL
                          	  return(obj)
                          }
                          )

setGeneric(name = "flAt", def = function(obj, idx) { standardGeneric("flAt") })

setMethod(f = "flAt", signature = "DensityData",
                          definition = function(obj, idx)
                          {
                              newData <- obj@data[idx, ]
                              rownames(newData) <- NULL
                          	  initialize(obj, data = newData)
                          }
                          )

setGeneric(name = "flSize", def = function(obj) { standardGeneric("flSize") })

setMethod(f = "flSize", signature = "DensityData",
                          definition = function(obj)
                          {
                          	  nrow(obj@data)
                          }
                          )

setGeneric(name = "flConcat", def = function(obj1, obj2) { standardGeneric("flConcat") })

setMethod(f = "flConcat", signature = "DensityData",
                          definition = function(obj1, obj2)
                          {
                              DensityData(data = rbind(obj1@data, obj2@data), numFeatures = obj1@numFeatures)
                          }
                          )

#' @export
flNormalizePopulationName <- function(population)
{
    populationNormalized <- stringr::str_replace_all(population, '\\+', 'p')
    populationNormalized <- stringr::str_replace_all(populationNormalized, '-', 'n')
    populationNormalized <- tolower(stringr::str_replace_all(populationNormalized, '[^a-zA-Z0-9]+', ''))
    populationNormalized
}

#' @export
flEstimateDensity <- function(data, n)
{
  dens <- density(data[which(!is.na(data))], n = n, from = 0)
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

#' @export
flDerivativeDtwDistanceMatrix <- function(densA, densB)
{
	lenA <- length(densA$y)
	lenB <- length(densB$y)

    if(lenA != lenB)
    {
        stop('Densities A and B must have same length.')
    }

    difA <- sapply(2:(lenA-1), function(i) (densA$y[i] - densA$y[i-1]) + ((densA$y[i+1] - densA$y[i-1]) / 2) / 2)
    # difA <- sapply(2:(lenA-1), function(i) (densA$y[i+1] - densA$y[i-1]) / (densA$x[i+1] - densA$x[i-1]))
    difA <- c(difA[1], difA, difA[lenA-2])

    difB <- sapply(2:(lenB-1), function(i) (densB$y[i] - densB$y[i-1]) + ((densB$y[i+1] - densB$y[i-1]) / 2) / 2)
    # difB <- sapply(2:(lenB-1), function(i) (densB$x[i+1] - densB$y[i-1]) / (densB$x[i+1] - densB$x[i-1]))
    difB <- c(difB[1], difB, difB[lenB-2])

    proxy::dist(difA, difB)
}

#' @export
flDtwMain <- function(densA, densB, ...)
{
    ddtw <- dtw::dtw(flDerivativeDtwDistanceMatrix(densA, densB), step.pattern = dtw::typeIds, ...)

    ddtw$query <- as.double(densA$y)
    ddtw$reference <- as.double(densB$y)

    ddtw
}

#' @export
flSelectPrototypes <- function(densdat, k)
{	
	D <- dist(flGetDensity(densdat)$y)

    protoIdx <- cluster::pam(D, k = k)$medoids

	return(protoIdx)
}

#' @export
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

#' @export
flPredictThresholds <- function(densdat, protoIdx)
{
	n <- flSize(densdat)
	nProto <- length(protoIdx)

	testIdx <- setdiff(1:n, protoIdx)
	nTest <- length(testIdx)

	D <- as.matrix(dist(flGetDensity(densdat)$y))

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

#' @export
flPlotDensThresh <- function(dens, thresh = NULL, predicted = NULL)
{
    p <- plot(
         dens$x,
         dens$y,
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

    return(p)
}