# Converts flowDensity gates from the OneStudy panel into evaluation files for flowLearn

library(stringr)
library(flowCore)
library(flowDensity)
library(flowLearn)

source("flowLearn.R")

setwd('onestudy')
source("OneStudy-Panel1.R")
setwd('..')

dir.create("trainingFiles.onestudy", showWarnings = FALSE)
dir.create("trainingFiles.onestudy/fcs", showWarnings = FALSE)

numFeatures <- 512

samples <- sampleNames(fs)

for (sample in samples)
{

    sampleTrainGates <- list()
    parentFrames <- list()

    # =============

    parentFrames[["CD45"]] <- fs.sngl[[sample]]
    sampleTrainGates[["CD45"]] <- TrainingGate(
                 parentName = "Singlets",
                 channelA = 12,
                 channelB = 3,
                 thresholdALow = CD45[[sample]]@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD45"]]) %in% CD45[[sample]]@index
                 )
    saveRDS(parentFrames[["CD45"]]@exprs[,c(12,3)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd45.parent.fcs.rds'))

    # =============

    parentFrames[["PBMC"]] <- fs.45[[sample]]
    sampleTrainGates[["PBMC"]] <- TrainingGate(
                 parentName = "CD45",
                 channelA = 7,
                 channelB = 3,
                 thresholdALow = NaN,
                 thresholdAHigh = PBMC[[sample]]@gates[1],
                 thresholdBLow = PBMC[[sample]]@gates[2],
                 thresholdBHigh = NaN,
                 negate = T,
                 gateAssignments = 1:nrow(parentFrames[["PBMC"]]) %in% PBMC[[sample]]@index
                 )
    saveRDS(parentFrames[["PBMC"]]@exprs[,c(7,3)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.pbmc.parent.fcs.rds'))

    # =============

    # Lymphocytes not available

    # =============

    parentFrames[["Bcell"]] <- fs.lymph[[sample]]
    sampleTrainGates[["Bcell"]] <- TrainingGate(
                 parentName = "Lymphocytes",
                 channelA = 6,
                 channelB = 10,
                 thresholdALow = Bcells[[sample]]$quad2@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = Bcells[[sample]]$quad2@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["Bcell"]]) %in% Bcells[[sample]]$quad2@index
                 )
    saveRDS(parentFrames[["Bcell"]]@exprs[,c(6,10)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.bcell.parent.fcs.rds'))

    # =============

    parentFrames[["Tcell"]] <- fs.lymph[[sample]]
    sampleTrainGates[["Tcell"]] <- TrainingGate(
                 parentName = "Lymphocytes",
                 channelA = 6,
                 channelB = 10,
                 thresholdALow = NaN,
                 thresholdAHigh = Bcells[[sample]]$quad3@gates[1],
                 thresholdBLow = Bcells[[sample]]$quad3@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["Tcell"]]) %in% Bcells[[sample]]$quad3@index
                 )
    saveRDS(parentFrames[["Tcell"]]@exprs[,c(6,10)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.tcell.parent.fcs.rds'))

    # =============

    parentFrames[["CD4"]] <- fs.tcell[[sample]]
    sampleTrainGates[["CD4"]] <- TrainingGate(
                 parentName = "Tcell",
                 channelA = 8,
                 channelB = 9,
                 thresholdALow = CD4.8[[sample]][['cd4']]@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = CD4.8[[sample]][['cd4']]@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD4"]]) %in% CD4.8[[sample]][['cd4']]@index
                 )
    saveRDS(parentFrames[["CD4"]]@exprs[,c(8,9)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd4.parent.fcs.rds'))

    # =============

    parentFrames[["CD8"]] <- fs.tcell[[sample]]
    sampleTrainGates[["CD8"]] <- TrainingGate(
                 parentName = "Tcell",
                 channelA = 8,
                 channelB = 9,
                 thresholdALow = NaN,
                 thresholdAHigh = CD4.8[[sample]][['cd8']]@gates[1],
                 thresholdBLow = CD4.8[[sample]][['cd8']]@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD8"]]) %in% CD4.8[[sample]][['cd8']]@index
                 )
    saveRDS(parentFrames[["CD8"]]@exprs[,c(8,9)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd8.parent.fcs.rds'))

    # =============

    parentFrames[["DN"]] <- fs.tcell[[sample]]
    sampleTrainGates[["DN"]] <- TrainingGate(
                 parentName = "Tcell",
                 channelA = 8,
                 channelB = 9,
                 thresholdALow = NaN,
                 thresholdAHigh = CD4.8[[sample]][['dn']]@gates[1],
                 thresholdBLow = NaN,
                 thresholdBHigh = CD4.8[[sample]][['dn']]@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["DN"]]) %in% CD4.8[[sample]][['dn']]@index
                 )
    saveRDS(parentFrames[["DN"]]@exprs[,c(8,9)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.dn.parent.fcs.rds'))

    # =============

    x <- rotate.data(fs.pbmc[[sample]], chans = c(4,7), theta = pi/14)$data
    parentFrames[["CD14"]] <- x
    sampleTrainGates[["CD14"]] <- TrainingGate(
                 parentName = "PBMC",
                 channelA = 4,
                 channelB = 7,
                 thresholdALow = NaN,
                 thresholdAHigh = NaN,
                 thresholdBLow = CD14[[sample]]@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD14"]]) %in% CD14[[sample]]@index
                 )
    saveRDS(parentFrames[["CD14"]]@exprs[,c(4,7)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd14.parent.fcs.rds'))

    # =============

    parentFrames[["CD14++CD16-"]] <- fs.14[[sample]]
    sampleTrainGates[["CD14++CD16-"]] <- TrainingGate(
                 parentName = "CD14",
                 channelA = 4,
                 channelB = 7,
                 thresholdALow = NaN,
                 thresholdAHigh = CD14.hi[[sample]]$quad3@gates[1],
                 thresholdBLow = CD14.hi[[sample]]$quad3@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD14++CD16-"]]) %in% CD14.hi[[sample]]$quad3@index
                 )
    saveRDS(parentFrames[["CD14++CD16-"]]@exprs[,c(4,7)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd14ppcd16n.parent.fcs.rds'))

    # =============

    parentFrames[["CD14++CD16+"]] <- fs.14[[sample]]
    sampleTrainGates[["CD14++CD16+"]] <- TrainingGate(
                 parentName = "CD14",
                 channelA = 4,
                 channelB = 7,
                 thresholdALow = CD14.hi[[sample]]$quad4@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = CD14.hi[[sample]]$quad4@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD14++CD16+"]]) %in% CD14.hi[[sample]]$quad4@index
                 )
    saveRDS(parentFrames[["CD14++CD16+"]]@exprs[,c(4,7)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd14ppcd16p.parent.fcs.rds'))

    # =============

    parentFrames[["CD14+CD16-"]] <- fs.14[[sample]]
    sampleTrainGates[["CD14+CD16-"]] <- TrainingGate(
                 parentName = "CD14",
                 channelA = 4,
                 channelB = 7,
                 thresholdALow = NaN,
                 thresholdAHigh = CD14.hi[[sample]]$quad1@gates[1],
                 thresholdBLow = NaN,
                 thresholdBHigh = CD14.hi[[sample]]$quad1@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD14+CD16-"]]) %in% CD14.hi[[sample]]$quad1@index
                 )
    saveRDS(parentFrames[["CD14+CD16-"]]@exprs[,c(4,7)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd14pcd16n.parent.fcs.rds'))

    # =============

    parentFrames[["CD14+CD16+"]] <- fs.14[[sample]]
    sampleTrainGates[["CD14+CD16+"]] <- TrainingGate(
                 parentName = "CD14",
                 channelA = 4,
                 channelB = 7,
                 thresholdALow = CD14.hi[[sample]]$quad2@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = CD14.hi[[sample]]$quad2@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD14+CD16+"]]) %in% CD14.hi[[sample]]$quad2@index
                 )
    saveRDS(parentFrames[["CD14+CD16+"]]@exprs[,c(4,7)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd14pcd16p.parent.fcs.rds'))

    # =============

    parentFrames[["NK"]] <- fs.lymph[[sample]]
    sampleTrainGates[["NK"]] <- TrainingGate(
                 parentName = "Lymphocytes",
                 channelA = 5,
                 channelB = 10,
                 thresholdALow = CD56[[sample]]$nk@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = CD56[[sample]]$nk@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["NK"]]) %in% CD56[[sample]]$nk@index
                 )
    saveRDS(parentFrames[["NK"]]@exprs[,c(5,10)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.nk.parent.fcs.rds'))

    # =============

    parentFrames[["NKT"]] <- fs.lymph[[sample]]
    sampleTrainGates[["NKT"]] <- TrainingGate(
                 parentName = "Lymphocytes",
                 channelA = 5,
                 channelB = 10,
                 thresholdALow = CD56[[sample]]$nkt@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = CD56[[sample]]$nkt@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["NKT"]]) %in% CD56[[sample]]$nkt@index
                 )
    saveRDS(parentFrames[["NKT"]]@exprs[,c(5,10)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.nkt.parent.fcs.rds'))

    # =============

    parentFrames[["CD56Bright"]] <- fs.nk[[sample]]
    sampleTrainGates[["CD56Bright"]] <- TrainingGate(
                 parentName = "NK",
                 channelA = 5,
                 channelB = 4,
                 thresholdALow = CD56.nk[[sample]]@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD56Bright"]]) %in% CD56.nk[[sample]]@index
                 )
    saveRDS(parentFrames[["CD56Bright"]]@exprs[,c(5,4)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd56bright.parent.fcs.rds'))

    # =============

    parentFrames[["CD64++CD16+"]] <- fs.14[[sample]]
    sampleTrainGates[["CD64++CD16+"]] <- TrainingGate(
                 parentName = "CD14",
                 channelA = 11,
                 channelB = 4,
                 thresholdALow = CD64.14[[sample]]@gates[1],
                 thresholdAHigh = NaN,
                 thresholdBLow = CD64.14[[sample]]@gates[2],
                 thresholdBHigh = NaN,
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD64++CD16+"]]) %in% CD64.14[[sample]]@index
                 )
    saveRDS(parentFrames[["CD64++CD16+"]]@exprs[,c(11,4)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd64ppcd16p.parent.fcs.rds'))

    # =============

    parentFrames[["CD56BrightCD16-"]] <- fs.56[[sample]]
    sampleTrainGates[["CD56BrightCD16-"]] <- TrainingGate(
                 parentName = "CD56Bright",
                 channelA = 5,
                 channelB = 4,
                 thresholdALow = NaN,
                 thresholdAHigh = NaN,
                 thresholdBLow = NaN,
                 thresholdBHigh = CD56.16[[sample]]@gates[2],
                 negate = F,
                 gateAssignments = 1:nrow(parentFrames[["CD56BrightCD16-"]]) %in% CD56.16[[sample]]@index
                 )
    saveRDS(parentFrames[["CD56BrightCD16-"]]@exprs[,c(5,4)], file = paste0('trainingFiles.onestudy/fcs/', sample, '.cd56brightcd16n.parent.fcs.rds'))

    for (pop in names(sampleTrainGates))
    {
        tg <- sampleTrainGates[[pop]]
        pf <- parentFrames[[pop]]


        exprsChanA <- pf@exprs[, tg@channelA]
        # negA <- which(exprsChanA < 0)
        # if(length(negA) > 0)
        # {
        #     exprsChanA <- exprsChanA[-negA]
        # }

        exprsChanB <- pf@exprs[, tg@channelB]
        # negB <- which(exprsChanB < 0)
        # if(length(negB) > 0)
        # {
        #     exprsChanB <- exprsChanB[-negB]
        # }


        sampleTrainGates[[pop]]@densitiesA[[as.character(numFeatures)]] <- estimateDensity(exprsChanA, numFeatures)
        sampleTrainGates[[pop]]@densitiesB[[as.character(numFeatures)]] <- estimateDensity(exprsChanB, numFeatures)
    }

    saveRDS(sampleTrainGates, file = paste('trainingFiles.onestudy/', sample, '.rds', sep = ''))

}

tr <- list()

tr[['cd45']] <- readTrainFiles('trainingFiles.onestudy', 'CD45', numFeatures)
tr[['bcell']] <- readTrainFiles('trainingFiles.onestudy', 'Bcell', numFeatures)
tr[['pbmc']] <- readTrainFiles('trainingFiles.onestudy', 'PBMC', numFeatures)
tr[['tcell']] <- readTrainFiles('trainingFiles.onestudy', 'Tcell', numFeatures)
tr[['cd4']] <- readTrainFiles('trainingFiles.onestudy', 'CD4', numFeatures)
tr[['cd8']] <- readTrainFiles('trainingFiles.onestudy', 'CD8', numFeatures)
tr[['dn']] <- readTrainFiles('trainingFiles.onestudy', 'DN', numFeatures)
tr[['cd14']] <- readTrainFiles('trainingFiles.onestudy', 'CD14', numFeatures)
tr[['cd14ppcd16n']] <- readTrainFiles('trainingFiles.onestudy', 'CD14++CD16-', numFeatures)
tr[['cd14ppcd16p']] <- readTrainFiles('trainingFiles.onestudy', 'CD14++CD16+', numFeatures)
tr[['cd14pcd16n']] <- readTrainFiles('trainingFiles.onestudy', 'CD14+CD16-', numFeatures)
tr[['cd14pcd16p']] <- readTrainFiles('trainingFiles.onestudy', 'CD14+CD16+', numFeatures)
tr[['nk']] <- readTrainFiles('trainingFiles.onestudy', 'NK', numFeatures)
tr[['nkt']] <- readTrainFiles('trainingFiles.onestudy', 'NKT', numFeatures)
tr[['cd56bright']] <- readTrainFiles('trainingFiles.onestudy', 'CD56Bright', numFeatures)
tr[['cd64ppcd16p']] <- readTrainFiles('trainingFiles.onestudy', 'CD64++CD16+', numFeatures)
tr[['cd56brightcd16n']] <- readTrainFiles('trainingFiles.onestudy', 'CD56BrightCD16-', numFeatures)

save(list = c('tr'), file = 'trainingFiles.onestudy/tr.RData')

# Calculate distance matrices for channel A and B

dA <- list()
dB <- list()

cl <- makeCluster(detectCores(), type = "FORK")
lapply(names(tr), function(x) 
{
    dists <- dtwDistanceMatrices(tr[[x]], cl)

    dA[[x]] <<- dists$dA
    dB[[x]] <<- dists$dB
})
stopCluster(cl)

save(list = c('tr', 'dA', 'dB'), file = 'trainingFiles.onestudy/tr.dtwdists.RData')