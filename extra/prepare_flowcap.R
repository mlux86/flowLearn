library(openCyto)
library(readxl)
library(flowLearn)
library(parallel)

printf <- function(...) invisible(print(sprintf(...)))

isRect <- function(gt)
{
    (class(gt)[1] == 'polygonGate') && (length(unique(gt@boundaries[,1])) == 2) && (length(unique(gt@boundaries[,2])) == 2)
}

boundariesToThresholds <- function(b)
{
    tAL <- min(b[, 1])
    tAH <- max(b[, 1])
    tBL <- min(b[, 2])
    tBH <- max(b[, 2])

    list(thresholdALow = tAL, thresholdAHigh = tAH, thresholdBLow = tBL, thresholdBHigh = tBH)
}

getParentIndices <- function(gs, popPath)
{
    which(getIndices(gs, getParent(gs, popPath))) %in% which(getIndices(gs, popPath))
}

flowcapPath <- '/mnt/data/immunespace/flowCAP-III-data/'
gslistPath <- paste0(flowcapPath, 'gating/gated_data/manual/')
fcsPath <- paste0(flowcapPath, 'fcs/SeraCare/')
p <- '/home/mlux/flowlearn.traindat/flowlearn.traindat.flowcap.'

sampleExcel <- read_excel(paste0(fcsPath, '000_files.xlsx'))

panels <- c('bcell', 'tcell', 'DC', 'treg')

cl <- parallel::makeCluster(parallel::detectCores(), type = "FORK", outfile = "")

# for (panel in panels)
parSapply(cl, panels, function(panel)
{

    dir.create(paste0(p, panel), showWarnings = FALSE)
    dir.create(paste0(p, panel, '/eval'), showWarnings = FALSE)

    gsl <- load_gslist(paste0(gslistPath, 'gslist-', panel))

    # generate gating ground truth

    densdat <- new('DensityData')

    filenames <- NULL
    samples <- NULL
    centers <- NULL
    replicates <- NULL

    for (sampleIdx in 1:length(gsl))
    {

        gs <- gsl[[sampleIdx]]
        nodes <- getNodes(gs)
        fcsFilename <- gs@name

        filenames <- c(filenames, fcsFilename)
        center <- sampleExcel[sampleExcel$Filename == fcsFilename, 'Center']
        centers <- c(centers, center)
        replicate <- sampleExcel[sampleExcel$Filename == fcsFilename, 'Replicate']
        replicates <- c(replicates, replicate)
        samples <- c(samples, sampleExcel[sampleExcel$Filename == fcsFilename, 'Sample'])

        gatingInfos <- list()
        for (nodeIdx in 2:length(nodes))
        {
            popPath <- nodes[nodeIdx]
            pops <- strsplit(popPath, '/')
            popName <- tail(pops[[1]], 1)

            parentPath <- getParent(gs, popPath)
            parentPops <- strsplit(parentPath, '/')
            parentName <- tail(parentPops[[1]], 1)

            npn <- flNormalizePopulationName(popName)

            gt <- getGate(gs, popPath)
            if(isRect(gt))
            {
                parentFrame <- getData(gs, parentPath)
                channels <- colnames(gt@boundaries)
                channelIndices <- c(which(colnames(parentFrame) == channels[1]), which(colnames(parentFrame) == channels[2]))
                thresholds <- boundariesToThresholds(gt@boundaries)

                if (is.null(gatingInfos[[npn]]))
                {
                    gatingInfos[[npn]] <- new("GatingInfo", population = npn, parent = flNormalizePopulationName(parentName), channels = channelIndices)
                }

                eA <- parentFrame@exprs[, channelIndices[1]]
                eB <- parentFrame@exprs[, channelIndices[2]]
                idx <- eA > 0 & eB > 0
                densA <- flEstimateDensity(eA[idx], 512)
                densB <- flEstimateDensity(eB[idx], 512)

                densdat <- flAdd(densdat, fcsFilename, npn, 1, densA$x, densA$y, thresholds$thresholdALow, thresholds$thresholdAHigh)
                densdat <- flAdd(densdat, fcsFilename, npn, 2, densB$x, densB$y, thresholds$thresholdBLow, thresholds$thresholdBHigh)

                neg <- npn == 'live'

                evaldat <- list(parentFrame = parentFrame@exprs[, channelIndices], indices = getParentIndices(gs, popPath), negate = neg)

                saveRDS(evaldat, file = paste0(p, panel, '/eval/', fcsFilename, '.', npn, '.eval.rds'))
            }
        }

    }

    sampleMeta <- data.frame('fcs' = filenames, 'sample' = unlist(samples), 'center' = unlist(centers), 'replicate' = unlist(replicates), stringsAsFactors = F)

    save(densdat, gatingInfos, sampleMeta, file = paste0(p, panel, '/train_data.RData'))

})

parallel::stopCluster(cl)
