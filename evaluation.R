f1Score <- function(sample, population, predictedThreshA, predictedThreshB)
# Calculates precision, recall and F1 score.
#
# Args:
#   sample: The name of the sample for loading the parent population and true gate information.
#   population: The name of the child population.
#   predictedThreshA: A two-element vector containing predicted lower and upper thresholds for density A (or NaN if not exists)
#   predictedThreshB: A two-element vector containing predicted lower and upper thresholds for density B (or NaN if not exists)
#
# Returns:
#   A 5-element vector c(precision, recall, f1, proportion, predicted proportion)
{
    tg <- readRDS(paste('trainingFiles/', sample , sep = ''))
    gateAssignments <- tg[[population]]@gateAssignments

    predictedGateAssignments <- gate(sample, population, predictedThreshA, predictedThreshB, tg[[population]]@negate)

    precision <- sum(gateAssignments & predictedGateAssignments) / sum(predictedGateAssignments)
    recall <- sum(gateAssignments & predictedGateAssignments) / sum(gateAssignments)
    f1 <- 2 * precision * recall / (precision + recall)
      
    c(precision, recall, f1, sum(gateAssignments), sum(predictedGateAssignments))
}

evaluatePerformance <- function(tr, population, testIdx, predictedThreshA, predictedThreshB)
# Evaluates performance for a full population.
#
# Args:
#   tr: A LearningSet object containing n samples.
#   population: The name of the child population.
#   testIdx: Indexes of samples to be evaluated.
#   predictedThreshA: A n-by-2 matrix containing predicted lower and upper thresholds for density A (or NaN if not exists)
#   predictedThreshB: A n-by-2 matrix containing predicted lower and upper thresholds for density B (or NaN if not exists)
#
# Returns: 
#   A list with keys meanPerf and medianPerf, each containing
#   5-element vectors with the mean/median precision, recall, f1 scores, proportions, and predicted proportions.
{
    samples <- tr@samples

    cl <- parallel::makeCluster(parallel::detectCores())
    parallel::clusterExport(cl, c("gate", "f1Score", "samples", "testIdx", "population", "predictedThreshA", "predictedThreshB"),  envir = environment())

    tryCatch({

        numTest <- length(testIdx)

        perf <- t(parallel::parSapply(cl, 1:numTest, function(i) 
        {
            f1Score(samples[[testIdx[i]]], population, predictedThreshA[i,], predictedThreshB[i,])
        }))

        p1 <- apply(perf, 2, function(x) mean(na.omit(x)))
        p2 <- apply(perf, 2, function(x) median(na.omit(x)))
        printf("\tprecision\trecall\t\tf1\n")
        printf("mean\t%.2f\t\t%.2f\t\t%.2f\n", p1[1], p1[2], p1[3])
        printf("median\t%.2f\t\t%.2f\t\t%.2f\n", p2[1], p2[2], p2[3])

        return(list(perf = perf, meanPerf = p1, medianPerf = p2))
    }, error = function(e) {
        print(e)
    }, finally = {
        parallel::stopCluster(cl) 
    })    
}

evaluatePopulation <- function(population, numTrain, seed = NaN, preloaded = NULL)
# Evaluates a full population, specifically loads it (if necessary), selects
# prototypes, predicts thresholds and evaluates precision, recall and F1.
#
# Args:
#   population: The name of the child population.
#   numTrain: Number of prototypes to use (i.e. 2 for 1 per channel)
#   seed: Seed for permuting samples after loading.
#   preloaded: A list of preloaded data to avoid repeated loading of data in some scenarios.
#              List keys are tr, dA and dB containing a loaded LearningSet object, 
#              and distance matrices for channel A and B, respectively.
#
# Returns: 
#   Performance as returned by evaluatePerformance()
{

    if(is.null(preloaded))
    {
        load('trainingFiles/tr.dtwdists.RData')
    } else
    {
        tr <- preloaded$tr
        dA <- preloaded$dA
        dB <- preloaded$dB
    }

    popId <- normalizePopulationName(population)

    tr <- tr[[popId]]
    dA <- dA[[popId]]
    dB <- dB[[popId]]

    n <- nrow(tr@densYchanA)

    # shuffle

    # if(!is.nan(seed))
    # {
    #     set.seed(seed)
    # }
    # perm <- sample(n)
    # dA <- dA[perm, perm]
    # dB <- dB[perm, perm]
    # tr <- permTrainFiles(tr, perm = perm)

    # predict

    selectedPrototypes <- selectPrototypes(tr, dA, dB, numTrain)

    predictedThresholds <- predictThresholds(tr, selectedPrototypes)

    evaluatePerformance(tr, population, selectedPrototypes$testIdx, predictedThresholds$threshA, predictedThresholds$threshB)
}

