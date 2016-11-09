#' @export
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

#' @export
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