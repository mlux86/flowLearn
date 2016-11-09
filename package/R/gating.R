#' @export
gate <- function(sample, population, predictedThreshA, predictedThreshB, negate)
# Gates a given population using provided thresholds
#
# Args:
#   sample: The name of the sample for loading the parent population and true gate information.
#   population: The name of the child population.
#   predictedThreshA: A two-element vector containing predicted lower and upper thresholds for density A (or NaN if not exists)
#   predictedThreshB: A two-element vector containing predicted lower and upper thresholds for density B (or NaN if not exists)
#
# Returns:
#   A logical vector with the number cells entries, defining whether a cell is in the gated population or not
{
    populationNormalized <- normalizePopulationName(population)

    fcs.name <- str_replace(sample, '\\.rds', '')
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