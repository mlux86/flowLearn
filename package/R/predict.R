#' Aligns a threshold from a reference to a test density.
#'
#' @param densX Density x values of the test density.
#' @param densY Density y values of the test density.
#' @param refDensX Density x values of the reference density.
#' @param refDensY Density y values of the reference density.
#' @param refThreshold Threshold of the reference density.
#'
#' @return The predicted threshold in the test density.
alignThreshold <- function(densX, densY, refDensX, refDensY, refThreshold)
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

#' Selects prototypes for a given set of samples.
#'
#' For now, prototypes are selected randomly.
#'
#' @param tr The LearningSet object.
#' @param dA The DTW distances for channel A.
#' @param dB The DTW distances for channel B.
#' @param numTrain The number of prototypes to use (i.e. 2 for 1 per channel)
#' @param seed A random seed.
#'
#' @return A list with the following keys:
#' \itemize{
#'  \item{"testIdx": Test indices not containing the selected prototype indices.}
#'  \item{"prototypesA": Prototype indices for channel A.}
#'  \item{"prototypesB": Prototype indices for channel B.}
#' }
#'
#' @export
selectPrototypes <- function(tr, dA, dB, numTrain, seed = NaN)
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

#' Selects prototypes for a given set of samples.
#'
#' For now, prototypes are selected randomly.
#'
#' @param tr The LearningSet object.
#' @param dA The DTW distances for channel A.
#' @param dB The DTW distances for channel B.
#'  @param trainIdxA Prototype indices for channel A.
#'  @param trainIdxB Prototype indices for channel B.
#'
#' @return A list with the following keys:
#' \itemize{
#'  \item{"testIdx": Test indices not containing the selected prototype indices.}
#'  \item{"prototypesA": Prototype indices for channel A.}
#'  \item{"prototypesB": Prototype indices for channel B.}
#' }
#'
#' @export
selectFixedPrototypes <- function(tr, dA, dB, trainIdxA, trainIdxB)
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

#' Predicts thresholds for one population.
#'
#' @param tr The LearningSet object.
#' @param selectedPrototypes Prototypes as selected by the selectPrototypes() or selectFixedPrototypes() functions.
#'
#' @return A list with the following keys:
#' \itemize{
#'  \item{"threshA": A n-by-2 matrix with predicted thresholds for channel A. NaN if threshold does not exist.}
#'  \item{"threshB": A n-by-2 matrix with predicted thresholds for channel B. NaN if threshold does not exist.}
#' }
#'
#' @export
predictThresholds <- function(tr, selectedPrototypes)
{
	numTest <- length(selectedPrototypes$testIdx)

    # print(numTest)
    # print(selectedPrototypes)

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
                # print(testIdx)
                # print(protoIdx)
		        alignThreshold(tr@densXchanA[testIdx,], tr@densYchanA[testIdx,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshA[1,2]))
		{
		    predictedThreshA[, 2] <- parallel::parSapply(cl, 1:numTest, function(i)
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesA[i]
                # print(testIdx)
                # print(protoIdx)
		        alignThreshold(tr@densXchanA[testIdx,], tr@densYchanA[testIdx,], tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx, 2])
		    })
		}
		if(!is.nan(tr@threshB[1,1]))
		{
		    predictedThreshB[, 1] <- parallel::parSapply(cl, 1:numTest, function(i)
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesB[i]
                # print(testIdx)
                # print(protoIdx)
		        alignThreshold(tr@densXchanB[testIdx,], tr@densYchanB[testIdx,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 1])
		    })
		}
		if(!is.nan(tr@threshB[1,2]))
		{
		    predictedThreshB[, 2] <- parallel::parSapply(cl, 1:numTest, function(i)
		    {
		    	testIdx <- selectedPrototypes$testIdx[i]
		        protoIdx <- selectedPrototypes$prototypesB[i]
                # print(testIdx)
                # print(protoIdx)
		        alignThreshold(tr@densXchanB[testIdx,], tr@densYchanB[testIdx,], tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx, 2])
		    })
		}

		return(list(threshA = predictedThreshA, threshB = predictedThreshB))

	}, error = function(e) {
		print(e)
	}, finally = {
		parallel::stopCluster(cl)
	})

}
