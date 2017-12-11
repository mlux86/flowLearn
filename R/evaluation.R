#' Given a DensityData object with existing thresholds, returns a logical vector indicating cell memberships.
#'
#' This method gates each cell and indicates whether it belongs to the given target population or not. For that, the method needs the original expression matrix.
#'
#' @param densdat The full DensityData object with density and densities and set/predicted thresholds.
#' @param exprs A n*d expression matrix of the parent population where n is the number of cells and d is the number of channels. The channel indices correspond to the channelIdx values in densdat.
#' @param fcs The name of the FCS file to consider. This is needed to filter the relevant entries in densdat.
#' @param population The name of the population to gate. This is needed to filter the relevant entries in densdat.
#' @param negate Whether the gated population should be negated/inverted, i.e. not the cells within the gate but out of it are considered.
#'
#' @return A vector of length n, indicating cell membership to the population of interest.
#'
#' @export
#'
#' @examples
#' f <- unique(flData(flSampleDensdat)$fcs)[[1]]
#' flGetGateAssignments(flSampleDensdat, flSampleBcellEvaluationData[[f]]$parentExprs, f, 'bcell')
#'
#'
flGetGateAssignments <- function(densdat, exprs, fcs, population, negate = FALSE)
{
    predictedGateAssignments <- matrix(TRUE, nrow(exprs), 1)

    nChan <- ncol(exprs)

    for (i in 1:nChan)
    {
        tmp <- flFind(densdat, fcs = fcs, population = population, channelIdx = i)

        g.l <- flGetGate(tmp)[1]
        g.h <- flGetGate(tmp)[2]

        if (!is.na(g.l))
        {
            predictedGateAssignments <- predictedGateAssignments & (exprs[,i] > g.l)
        }
        if (!is.na(g.h))
        {
            predictedGateAssignments <- predictedGateAssignments & (exprs[,i] < g.h)
        }
    }

    if (negate)
    {
        predictedGateAssignments <- !predictedGateAssignments
    }

    predictedGateAssignments
}

#' Given a DensityData object with existing thresholds, calculates precision, recall, F1-score, true and predicted cell proportions for one FCS file.
#'
#' This method gates each cell and compares the result to a true vector of true cell assignments. It returns its performance in the form of various parameters.
#'
#' @param densdat The full DensityData object with density and densities and set/predicted thresholds.
#' @param fcs The name of the FCS file to consider. This is needed to filter the relevant entries in densdat.
#' @param population The name of the population to gate. This is needed to filter the relevant entries in densdat.
#' @param trueAssignments A vector of length <number of cells>, indicating true cell memberships of the target population.
#' @param parentExprs A n*d expression matrix of the parent population where n is the number of cells and d is the number of channels. The channel indices correspond to the channelIdx values in densdat.
#' @param negate Whether the gated population should be negated/inverted, i.e. not the cells within the gate but out of it are considered.
#'
#' @return A data.frame with precision, recall, f1, trueProportion, predictedProportion.
#'
#' @export
#'
#' @examples
#' # In this example, we evaluate true against the true gate assignments, hence it should give perfect performance, i.e. F_1 = 1
#'
#' f <- unique(flData(flSampleDensdat)$fcs)[[1]]
#' flEvalF1ScoreFCS(flSampleDensdat, f, 'bcell', flSampleBcellEvaluationData[[f]]$gateAssignments, flSampleBcellEvaluationData[[f]]$parentExprs, FALSE)
#'
#'
flEvalF1ScoreFCS <- function(densdat, fcs, population, trueAssignments, parentExprs, negate = FALSE)
{
    predictedGateAssignments <- flGetGateAssignments(densdat, parentExprs, fcs, population, negate)

    precision <- sum(trueAssignments & predictedGateAssignments) / sum(predictedGateAssignments)
    recall <- sum(trueAssignments & predictedGateAssignments) / sum(trueAssignments)

    if (precision + recall == 0)
    {
        f1 <- 0
    } else
    {
        f1 <- 2 * precision * recall / (precision + recall)
    }

    data.frame(precision = precision, recall = recall, f1 = f1, trueProportion = sum(trueAssignments) / length(trueAssignments), predictedProportion = sum(predictedGateAssignments) / length(predictedGateAssignments))
}

#' Identify samples outliers by the means of F1-scores.
#'
#' This method returns indices of samples with low F1-scores, using an outlier heuristic that considers all samples as outliers, when they are below the 25-percent-quantile - 1.5 * inter-quartile-range.
#'
#' @param f1s The list of F1-scores, one entry per sample.
#'
#' @return A vector of indices with low F1-scores.
#'
#' @importFrom stats quantile
#' @importFrom stats IQR
#' @export
#'
#' @examples
#' f1s <- runif(99, min = 0.42, max = 0.84)
#' f1s[100] <- 0
#' flIdentifyOutliersF1(f1s)
#'
flIdentifyOutliersF1 <- function(f1s)
{
    quantiles <- stats::quantile(f1s, na.rm = TRUE)
    minn <- quantiles[[2]] - stats::IQR(f1s, na.rm = TRUE) * 1.5
    which(f1s < minn)
}

