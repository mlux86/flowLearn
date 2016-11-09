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
    dtwObj <- dtwMain(densY, refDensY)

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
#       prototypesA: Prototype indices for channel A.
#       prototypesB: Prototype indices for channel B.
{

	if(numTrain %% 2 != 0)
	{
		stop('Number of training samples must be one or divisible by two')
	}

	n <- nrow(tr@densYchanA)

	K <- numTrain / 2 # for A and B

	# random training indices

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

	# test indices don't contain training indices

	trainIdx <- union(trainIdxA, trainIdxB)
	testIdx <- setdiff(1:n, trainIdx)

	# find nearest prototype assignments for each sample

	if(K > 1)
	{
	    labelsA <- t(apply(dA[,trainIdxA], 1, order))[, 1]
	    labelsB <- t(apply(dB[,trainIdxB], 1, order))[, 1]
	} else 
	{
	    labelsA <- matrix(1, n, 1)
	    labelsB <- matrix(1, n, 1)
	}

	prototypesA <- trainIdxA[labelsA]
	prototypesB <- trainIdxB[labelsB]

	return(list(testIdx = testIdx, prototypesA = prototypesA[testIdx], prototypesB = prototypesB[testIdx]))

}

# #' @export
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
#       prototypesA: Prototype indices for channel A.
#       prototypesB: Prototype indices for channel B.
{
	n <- nrow(tr@densYchanA)

	K <- length(trainIdxA)
	if(length(trainIdxB) != K)
	{
		stop('Number of prototypes in both channels must be the same!')
	}

	# test indices don't contain training indices

	trainIdx <- union(trainIdxA, trainIdxB)
	testIdx <- setdiff(1:n, trainIdx)

	# find nearest prototype assignments for each sample

	if(K > 1)
	{
	    labelsA <- t(apply(dA[,trainIdxA], 1, order))[, 1]
	    labelsB <- t(apply(dB[,trainIdxB], 1, order))[, 1]
	} else 
	{
	    labelsA <- matrix(1, n, 1)
	    labelsB <- matrix(1, n, 1)
	}

	prototypesA <- trainIdxA[labelsA]
	prototypesB <- trainIdxB[labelsB]

	return(list(testIdx = testIdx, prototypesA = prototypesA[testIdx], prototypesB = prototypesB[testIdx]))
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

	cl <- parallel::makeCluster(parallel::detectCores())
	parallel::clusterExport(cl, c("tr", "selectedPrototypes", "numTest"),  envir = environment())

	tryCatch({

		if(!is.nan(tr@threshA[1,1]))
		{
		    predictedThreshA[, 1] <- parallel::parSapply(cl, 1:numTest, function(i) 
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesA[i]
		        alignThreshold(tr@densXchanA[i,], tr@densYchanA[i,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshA[1,2]))
		{
		    predictedThreshA[, 2] <- parallel::parSapply(cl, 1:numTest, function(i) 
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesA[i]
		        print(protoIdx)
		        alignThreshold(tr@densXchanA[i,], tr@densYchanA[i,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 2])
		    })
		}
		if(!is.nan(tr@threshB[1,1]))
		{
		    predictedThreshB[, 1] <- parallel::parSapply(cl, 1:numTest, function(i) 
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesB[i]
		        print(protoIdx)
		        alignThreshold(tr@densXchanB[i,], tr@densYchanB[i,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshB[1,2]))
		{
		    predictedThreshB[, 2] <- parallel::parSapply(cl, 1:numTest, function(i) 
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesB[i]
		        print(protoIdx)
		        alignThreshold(tr@densXchanB[i,], tr@densYchanB[i,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 2])
		    })
		}

		return(list(threshA = predictedThreshA, threshB = predictedThreshB))

	}, error = function(e) {
		print(e)
	}, finally = {
		parallel::stopCluster(cl)	
	})

}
