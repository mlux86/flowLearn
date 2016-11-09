#' Gates a given population using provided thresholds
#'
#' @param sample The name of the sample for loading the parent population and true gate information.
#' @param population The name of the child population.
#' @param predictedThreshA A two-element vector containing predicted lower and upper thresholds for density A (or NaN if not exists)
#' @param predictedThreshB A two-element vector containing predicted lower and upper thresholds for density B (or NaN if not exists)
#'
#' @return A logical vector with the number cells entries, defining whether a cell is in the gated population or not
#'
#' @export
gate <- function(sample, population, predictedThreshA, predictedThreshB, negate)
{
    populationNormalized <- normalizePopulationName(population)

    fcs.name <- stringr::str_replace(sample, '\\.rds', '')
    exprs <- readRDS(paste('trainingFiles/fcs/', fcs.name, '.', populationNormalized, '.parent.fcs.rds', sep = ''))

    predictedGateAssignments <- matrix(T, nrow(exprs), 1)
    if (!is.nan(predictedThreshA[1]))
    {
        predictedGateAssignments <- predictedGateAssignments & exprs[,1] > predictedThreshA[1]
    }
    if (!is.nan(predictedThreshA[2]))
    {
        predictedGateAssignments <- predictedGateAssignments & exprs[,1] < predictedThreshA[2]
    }    
    if (!is.nan(predictedThreshB[1]))
    {
        predictedGateAssignments <- predictedGateAssignments & exprs[,2] > predictedThreshB[1]
    }
    if (!is.nan(predictedThreshB[2]))
    {
        predictedGateAssignments <- predictedGateAssignments & exprs[,2] < predictedThreshB[2]
    }

    if (negate)
    {
        predictedGateAssignments <- !predictedGateAssignments
    }

    predictedGateAssignments
}

#' Calculates precision, recall and F1 score.
#'
#' @param sample The name of the sample for loading the parent population and true gate information.
#' @param population The name of the child population.
#' @param predictedThreshA A two-element vector containing predicted lower and upper thresholds for density A (or NaN if not exists)
#' @param predictedThreshB A two-element vector containing predicted lower and upper thresholds for density B (or NaN if not exists)
#'
#' @return A 5-element vector c(precision, recall, f1, proportion, predicted proportion)
f1Score <- function(sample, population, predictedThreshA, predictedThreshB)
{
    tg <- readRDS(paste('trainingFiles/', sample , sep = ''))
    gateAssignments <- tg[[population]]@gateAssignments

    predictedGateAssignments <- gate(sample, population, predictedThreshA, predictedThreshB, tg[[population]]@negate)

    precision <- sum(gateAssignments & predictedGateAssignments) / sum(predictedGateAssignments)
    recall <- sum(gateAssignments & predictedGateAssignments) / sum(gateAssignments)
    f1 <- 2 * precision * recall / (precision + recall)
      
    c(precision, recall, f1, sum(gateAssignments), sum(predictedGateAssignments))
}

#' Evaluates performance for a full population.
#'
#' @param tr A LearningSet object containing n samples.
#' @param population The name of the child population.
#' @param testIdx Indexes of samples to be evaluated.
#' @param predictedThreshA A n-by-2 matrix containing predicted lower and upper thresholds for density A (or NaN if not exists)
#' @param predictedThreshB A n-by-2 matrix containing predicted lower and upper thresholds for density B (or NaN if not exists)
#'
#' @return A list with keys meanPerf and medianPerf, each containing 5-element vectors with the 
#'         mean/median precision, recall, f1 scores, proportions, and predicted proportions.
evaluatePerformance <- function(tr, population, testIdx, predictedThreshA, predictedThreshB)
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

#' Evaluates a full population. 
#'
#' Specifically loads it (if necessary), selects prototypes, predicts thresholds and evaluates precision, recall and F1.
#'
#' @param population The name of the child population.
#' @param numTrain Number of prototypes to use (i.e. 2 for 1 per channel)
#' @param seed Seed for permuting samples after loading.
#' @param preloaded A list of preloaded data to avoid repeated loading of data in some scenarios.
#'                  List keys are tr, dA and dB containing a loaded LearningSet object, 
#'                  and distance matrices for channel A and B, respectively.
#'
#' @return Performance as returned by evaluatePerformance()
evaluatePopulation <- function(population, numTrain, seed = NaN, preloaded = NULL)
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

#' Reads proportions for a given population and return the mean proportion.
#'
#' @param tr The LearningSet object to calculate the statistic for.
#' @param population The population name.
#'
#' @return The mean proportion.
#'
#' @export
getMeanProportion <- function(tr, population)
{
    populationNormalized <- normalizePopulationName(population)

    n <- length(tr[[populationNormalized]]@samples)

    cl <- parallel::makeCluster(parallel::detectCores())
    parallel::clusterExport(cl, c("tr", "population", "populationNormalized"),  envir = environment())

    tryCatch({
        props <- t(parallel::parSapply(cl, 1:n, function(i) 
        {
            tg <- readRDS(paste0('trainingFiles/', tr[[populationNormalized]]@samples[[i]]))

            gateAssignments <- tg[[population]]@gateAssignments     

            sum(gateAssignments) / length(gateAssignments)
        }))

        return(mean(props))
    }, error = function(e) {
        print(e)
    }, finally = {
        parallel::stopCluster(cl) 
    }) 
}