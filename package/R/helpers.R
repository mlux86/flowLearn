#' Normalizes a population name.
#'
#' First, all '+' and '-' signs are converted to 'p' and 'n', respectively.
#' Second, All non-characters and non-numbers are removed.
#' Third, everything is converted to lower case.
#'
#' @param population The unnormalized population name string.
#'
#' @return The normalized population string.
#'
#' @export
normalizePopulationName <- function(population)
{
    populationNormalized <- stringr::str_replace_all(population, '\\+', 'p')
    populationNormalized <- stringr::str_replace_all(populationNormalized, '-', 'n')
    populationNormalized <- tolower(stringr::str_replace_all(populationNormalized, '[^a-zA-Z0-9]+', ''))
    populationNormalized
}

#' Estimates a smoothed density from a given vector of numbers.
#'
#' @param data The numbers to calculated the density for.
#' @param n Number of features for the density.
#'
#' @return The estimated density.
#'
#' @export
estimateDensity <- function(data, n)
{
  dens <- density(data[which(!is.na(data))], n = n)
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  return(dens)
}

#' Rotates FCS data by theta.
#'
#' @param data The flowFrame to rotate or the ecpression matrix.
#' @param chans Only rotate the given channels.
#' @param theta Amount of rotation to apply.
#'
#' @return A rotated version of the input.
#'
#' @export
rotateData <- function(data, chans=NULL, theta=NULL)
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
