# See http://bioinfosrv1.bccrc.ca/index.php/FlowLearn for explicit documentation

TrainingGate <- setClass("TrainingGate",
                         slots = c(
                                   parentName = "character",
                                   channelA = "numeric",
                                   channelB = "numeric",
                                   thresholdALow = "numeric",
                                   thresholdAHigh = "numeric",
                                   thresholdBLow = "numeric",
                                   thresholdBHigh = "numeric",
                                   negate = "logical",
                                   gateAssignments = "logical",
                                   densitiesA = "list",
                                   densitiesB = "list"
                                   )
                         )

LearningSet <- setClass("LearningSet",
                         slots = c(
                                   densXchanA = "matrix",
                                   densXchanB = "matrix",
                                   densYchanA = "matrix",
                                   densYchanB = "matrix",
                                   threshA = "matrix",
                                   threshB = "matrix",
                                   samples = "list"
                                   )
                         )

readTrainFiles <- function(folder, population, numFeatures = 512)
# Loads training data for one population.
#
# Args:
#   folder: Folder containing the data.
#   population: The population to load.
#   numFeatures: Number of features to use.
#
# Returns:
#   A LearningSet object that contains the requested data.
{
    files <- list.files(path = folder, full.names = T, recursive = F, pattern = '*.rds')
    filesNames <- list.files(path = folder, full.names = F, recursive = F, pattern = '*.rds')

    densXchanA <- matrix(NaN, ncol = numFeatures, nrow = length(files))
    densYchanA <- matrix(NaN, ncol = numFeatures, nrow = length(files))
    densXchanB <- matrix(NaN, ncol = numFeatures, nrow = length(files))
    densYchanB <- matrix(NaN, ncol = numFeatures, nrow = length(files))
    threshA <- matrix(NaN, ncol = 2, nrow = length(files))
    threshB <- matrix(NaN, ncol = 2, nrow = length(files))

    samples <- list()

    i <- 1
    for (f in files)
    {
        s <- readRDS(f)
        samples[[i]] <- filesNames[[i]]
        densXchanA[i, ] <- s[[population]]@densitiesA[[as.character(numFeatures)]]$x
        densYchanA[i, ] <- s[[population]]@densitiesA[[as.character(numFeatures)]]$y
        densXchanB[i, ] <- s[[population]]@densitiesB[[as.character(numFeatures)]]$x
        densYchanB[i, ] <- s[[population]]@densitiesB[[as.character(numFeatures)]]$y
        threshA[i, ] <- c(s[[population]]@thresholdALow, s[[population]]@thresholdAHigh)
        threshB[i, ] <- c(s[[population]]@thresholdBLow, s[[population]]@thresholdBHigh)

        i <- i + 1
    }

    LearningSet(        
        densXchanA = densXchanA,
        densYchanA = densYchanA,
        densXchanB = densXchanB,
        densYchanB = densYchanB,
        threshA = threshA,
        threshB = threshB,
        samples = samples
    )
}

permTrainFiles <- function(tr, perm = NaN, seed = NaN, subsampleRatio = NaN)
# Permutes a LearningSet object (as loaded by readTrainFiles).
#
# Args:
#   tr: The LearningSet object.
#   perm: The permutation to use (use the sample function).
#   seed: Random seed to use (ignored when perm is given).
#   subsampleRatio: Only return this proportion of the given data (i.e. for bootstrapping).
#                   (ignored if perm is given)
#
# Returns:
#   A LearningSet object that was permuted by the given parameters.
{

    n <- length(tr@samples)
    
    if (is.nan(perm[[1]]))
    {
        if (!is.nan(seed))
        {
            set.seed(seed)
        }

        perm <- sample(n)

        if (!is.nan(subsampleRatio))
        {
            perm <- perm[1:round(subsampleRatio * n)]
        }
    }

    tr@densXchanA <- tr@densXchanA[perm, ]
    tr@densYchanA <- tr@densYchanA[perm, ]
    tr@densXchanB <- tr@densXchanB[perm, ]
    tr@densYchanB <- tr@densYchanB[perm, ]
    tr@threshA <- tr@threshA[perm, ]
    tr@threshB <- tr@threshB[perm, ]
    tr@samples <- tr@samples[perm]

    tr
}

dtwDistanceMatrices <- function(tr, cl)
# Computes DTW distance matrices for channel A and B, respectively.
# Each matrix is of size nxn where n is the number of samples.
# This implementation is parallelized. Only the upper triangle matrix is computed and then mirrored.
#
# Args:
#   tr: A LearningSet object.
#   cl: A compute cluster as obtained by parallel::makeCluster()
#
# Returns:
#   A list with keys dA and dB indicating the distance matrices for channel A and B.
{

    n <- length(tr@samples)

    # channel A

    m <- parSapply(cl, 2:n, # start at 2 for upper triangle only
        function(i.1) 
        {
            sapply(1:(i.1-1), function(i.2) myDtw(tr@densYchanA[i.1,], tr@densYchanA[i.2,], distance.only = T)$distance)
        })

    dA <- matrix(0, n, n)
    dA[lower.tri(dA, diag=FALSE)] <- unlist(m)
    dA <- dA + t(dA)

    # channel B

    m <- parSapply(cl, 2:n, # start at 2 for upper triangle only
        function(i.1) 
        {
            sapply(1:(i.1-1), function(i.2) myDtw(tr@densYchanB[i.1,], tr@densYchanB[i.2,], distance.only = T)$distance)
        })

    dB <- matrix(0, n, n)
    dB[lower.tri(dB, diag=FALSE)] <- unlist(m)
    dB <- dB + t(dB)

    list(dA = dA, dB = dB)
}