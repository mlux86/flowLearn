# A helper function for the C-style printf() function.
printf <- function(...) invisible(cat(sprintf(...)))

#' @export
normalizePopulationName <- function(population)
# Normalizes a population name.
# First, all '+' and '-' signs are converted to 'p' and 'n', respectively.
# Second, All non-characters and non-numbers are removed.
# Third, everything is converted to lower case.
#
# Args:
#   population: The unnormalized population name string.
#
# Returns:
#   The normalized population string.
{
    populationNormalized <- stringr::str_replace_all(population, '\\+', 'p')
    populationNormalized <- stringr::str_replace_all(populationNormalized, '-', 'n')
    populationNormalized <- tolower(stringr::str_replace_all(populationNormalized, '[^a-zA-Z0-9]+', ''))
    populationNormalized
}

#' @export
estimateDensity <- function(data, n)
# Estimates a smoothed density from a given vector of numbers.
#
# Args:
#   data: The numbers to calculated the density for.
#   n: Number of features for the density.
#
# Returns:
#   The estimated density.
{
  dens <- density(data[which(!is.na(data))], n = n)
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

#' @export
rotateData <- function(data, chans=NULL, theta=NULL)
# Rotates FCS data by theta.
#
# Args:
#   data: The flowFrame to rotate or the ecpression matrix.
#   chans: Only rotate the given channels.
#   theta: Amount of rotation to apply.
#
# Returns:
#   A rotated version of the input.
{
    if (class(data)== "flowFrame" & !is.null(chans))
    {
        dataNew <- exprs(data)[,chans]
        if (is.null(theta))
        {
            regSlope <- atan(lm(dataNew[,1] ~ dataNew[,2])$coefficients[2])
            theta <- pi/2 - regSlope
        }
        dataNew <- dataNew %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
        exprs(data)[,chans] <- dataNew
    }else{
        data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    }
    return(list(data=data,theta=theta))
}

#' @export
getMeanProportion <- function(tr, population)
# Reads proportions for a given population and return the mean proportion.
#
# Args:
#   tr: The LearningSet object to calculate the statistic for.
#   population: The population name.
#
# Returns:
#   The mean proportion.
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