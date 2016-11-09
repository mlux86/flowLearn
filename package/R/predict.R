alignThreshold <- function(densX, densY, refDensX, refDensY, refThreshold)
# Aligns a threshold from a reference to a test density.
#
# Args:
#   densX: Density x values of the test density.
#   densY: Density y values of the test density.
#   refDensX: Density x values of the reference density.
#   refDensY: Density y values of the reference density.
#   refThreshold: Threshold of the reference density.
#
# Returns:
#   The predicted threshold in the test density.
{
    dtwObj <- myDtw(densY, refDensY)

    # find x index nearest to threshold in reference

    nearestRefXIdx <- which.min(abs( refDensX - refThreshold ))    

    # identify mapped indexes on aligned density

    idx <- dtwObj$index1[dtwObj$index2 == nearestRefXIdx]

    if (length(idx) > 1)
    {
        thresh <- mean(densX[idx])
    } else
    {
        thresh <- densX[idx]
    }

    thresh
}

#' @export
selectPrototypes <- function(tr, dA, dB, numTrain, seed = NaN)
# Selects prototypes for a given set of samples.
# For now, prototypes are selected randomly.
#
# Args:
#   tr: The LearningSet object.
#   dA: The DTW distances for channel A.
#   dB: The DTW distances for channel B.
#   numTrain: The number of prototypes to use (i.e. 2 for 1 per channel)
#   seed: A random seed.
#
# Returns:
#   A list with the following keys:
#   	testIdx: Test indices not containing the selected prototype indices.
#       trainIdxA: Prototype indices for channel A.
#       trainIdxB: Prototype indices for channel B.
#       labelsA: A vector of length equal to the number of samples and containg the 
#                prototype index for each test sample (channel A).
#       labelsB: A vector of length equal to the number of samples and containg the 
#                prototype index for each test sample (channel B).
{

	if(numTrain %% 2 != 0)
	{
		stop('Number of training samples must be one or divisible by two')
	}

	n <- nrow(tr@densYchanA)

	K <- numTrain / 2 # for A and B

    if(!is.nan(seed))
    {
        set.seed(seed)
    }
	trainIdxA <- sample(n, K)
    if(!is.nan(seed))
    {
        set.seed(seed + 1)
    }	
	trainIdxB <- sample(n, K)

	if(K > 1)
	{
	    labelsA <- t(apply(dA[,trainIdxA], 1, order))[, 1]
	    labelsB <- t(apply(dB[,trainIdxB], 1, order))[, 1]
	} else 
	{
	    labelsA <- matrix(1, n, 1)
	    labelsB <- matrix(1, n, 1)
	}

	testIdxA <- (1:n)[-trainIdxA]
	testIdxB <- (1:n)[-trainIdxB]

	trainIdx <- union(trainIdxA, trainIdxB)
	testIdx <- setdiff(1:n, trainIdx)

	return(list(testIdx = testIdx, trainIdxA = trainIdxA, trainIdxB = trainIdxB, labelsA = labelsA, labelsB = labelsB))

}

#' @export
selectFixedPrototypes <- function(tr, dA, dB, trainIdxA, trainIdxB)
# Selects prototypes for a given set of samples.
# Fixes the prototype indices and only calculates test indices to be used.
#
# Args:
#   tr: The LearningSet object.
#   dA: The DTW distances for channel A.
#   dB: The DTW distances for channel B.
#   trainIdxA: Prototype indices for channel A.
#   trainIdxB: Prototype indices for channel B.
#
# Returns:
#   A list with the following keys:
#   	testIdx: Test indices not containing the selected prototype indices.
#       trainIdxA: Prototype indices for channel A.
#       trainIdxB: Prototype indices for channel B.
#       labelsA: A vector of length equal to the number of samples and containg the 
#                prototype index for each test sample (channel A).
#       labelsB: A vector of length equal to the number of samples and containg the 
#                prototype index for each test sample (channel B).
{
	n <- nrow(tr@densYchanA)

	K <- length(trainIdxA)
	if(length(trainIdxB) != K)
	{
		stop('Number of prototypes in both channels must be the same!')
	}

	if(K > 1)
	{
	    labelsA <- t(apply(dA[,trainIdxA], 1, order))[, 1]
	    labelsB <- t(apply(dB[,trainIdxB], 1, order))[, 1]
	} else 
	{
	    labelsA <- matrix(1, n, 1)
	    labelsB <- matrix(1, n, 1)
	}	

	testIdxA <- (1:n)[-trainIdxA]
	testIdxB <- (1:n)[-trainIdxB]

	trainIdx <- union(trainIdxA, trainIdxB)
	testIdx <- setdiff(1:n, trainIdx)

	return(list(testIdx = testIdx, trainIdxA = trainIdxA, trainIdxB = trainIdxB, labelsA = labelsA, labelsB = labelsB))
}

#' @export
predictThresholds <- function(tr, selectedPrototypes)
# Predicts thresholds for one population.
#
# Args:
#   tr: The LearningSet object.
#   selectedPrototypes: Prototypes as selected by the selectPrototypes() or selectFixedPrototypes() functions.
#
# Returns:
#   A list with the following keys:
#   	threshA: A n-by-2 matrix with predicted thresholds for channel A. NaN if threshold does not exist.
#       threshB: A n-by-2 matrix with predicted thresholds for channel B. NaN if threshold does not exist.
{
	numTest <- length(selectedPrototypes$testIdx)

	predictedThreshA <- matrix(NaN, numTest, 2)
	predictedThreshB <- matrix(NaN, numTest, 2)

	cl <- makeCluster(detectCores(), type = "FORK")

	tryCatch({

		if(!is.nan(tr@threshA[1,1]))
		{
		    predictedThreshA[, 1] <- parSapply(cl, selectedPrototypes$testIdx, function(i) 
		    {
		        protoIdx <- selectedPrototypes$trainIdxA[selectedPrototypes$labelsA[i]]
		        alignThreshold(tr@densXchanA[i,], tr@densYchanA[i,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshA[1,2]))
		{
		    predictedThreshA[, 2] <- parSapply(cl, selectedPrototypes$testIdx, function(i) 
		    {
		        protoIdx <- selectedPrototypes$trainIdxA[selectedPrototypes$labelsA[i]]
		        alignThreshold(tr@densXchanA[i,], tr@densYchanA[i,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 2])
		    })
		}
		if(!is.nan(tr@threshB[1,1]))
		{
		    predictedThreshB[, 1] <- parSapply(cl, selectedPrototypes$testIdx, function(i) 
		    {
		        protoIdx <- selectedPrototypes$trainIdxB[selectedPrototypes$labelsB[i]]
		        alignThreshold(tr@densXchanB[i,], tr@densYchanB[i,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshB[1,2]))
		{
		    predictedThreshB[, 2] <- parSapply(cl, selectedPrototypes$testIdx, function(i) 
		    {
		        protoIdx <- selectedPrototypes$trainIdxB[selectedPrototypes$labelsB[i]]
		        alignThreshold(tr@densXchanB[i,], tr@densYchanB[i,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 2])
		    })
		}

		return(list(threshA = predictedThreshA, threshB = predictedThreshB))

	}, error = function(e) {
		print(e)
	}, finally = {
		stopCluster(cl)	
	})

}
