# Converts flowDensity gates from the IMPC BM panel into evaluation files for flowLearn

library(stringr)
library(flowCore)
library(flowDensity)

source("flowLearn.R")

dir.create("trainingFiles.impc", showWarnings = FALSE)
dir.create("trainingFiles.impc/fcs", showWarnings = FALSE)

numFeatures <- 512

cleanPath <- paste("/mnt/f/Brinkman group/current/Albina/Markus/After_Clean")
thresholdPath <- paste("/mnt/f/Brinkman group/current/Albina/Markus/Gating_Thresholds_Updated")

samples <- list.files(path = cleanPath, full.names = T, recursive = T, pattern = '*.Rdata')

cl <- makeCluster(detectCores(), type = "FORK", outfile = "")

lst <- 1:length(samples)
# lst <- sample(length(samples), 100)

parSapply(cl, lst, function(idx) {

    s <- samples[[idx]]

    tryCatch({

        load(s)

        fpath <- f@description$FILENAME
        fname <- basename(fpath)

        sample <- str_replace(fname, "\\.fcs", "")
        genotype <- strsplit(fpath, split = "/")[[1]][[7]]

        print(paste(genotype, sample, sep = "/"))

        load(file = paste(thresholdPath, "/", genotype, "/Gthres_", fname, ".Rdata", sep = ""))

        sampleTrainGates <- list()
        parentFrames <- list()

        # Gating All Events to obtain the Singlets. Plotting SSC-A_FSC-W
        singlets.flowD.h <- flowDensity(f, channels = c(4, 3), position = c(NA, F), gates = c(NA, gthres[1]))
        singlets.flowD.l <- flowDensity(singlets.flowD.h, channels = c(4, 3), position = c(NA, T), gates = c(NA, gthres[2]))
        singlets.flowD.l@proportion <- (singlets.flowD.l@cell.count/nrow(f))*100
        singlets.flowD.l.ind <- singlets.flowD.l@index
        singlets <- getflowFrame(singlets.flowD.l)

        parentFrames[["Singlets"]] <- f
        sampleTrainGates[["Singlets"]] <- TrainingGate(
                     parentName = "",
                     channelA = 4,
                     channelB = 3,
                     thresholdALow = NaN,
                     thresholdAHigh = NaN,
                     thresholdBLow = gthres[2],
                     thresholdBHigh = gthres[1],
                     negate = F,
                     gateAssignments = 1:nrow(f) %in% singlets.flowD.l@index
                     )

        saveRDS(f@exprs[,c(4,3)], file = paste('trainingFiles.impc/fcs/', sample, '.singlets.parent.fcs.rds', sep = ''))

        # Gating Singlets to obtain the Live population. Plotting Live/Dead_SSC-A -----
        live.flowD <- flowDensity(singlets, channels = c(11,4), position = c(F,NA), gates = c(gthres[3], NA))
        live <- getflowFrame(live.flowD)

        parentFrames[["Live"]] <- singlets
        sampleTrainGates[["Live"]] <- TrainingGate(
                     parentName = "Singlets",
                     channelA = 11,
                     channelB = 4,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[2],
                     thresholdBLow = NaN,
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(singlets) %in% live.flowD@index
                     )

        saveRDS(singlets@exprs[,c(11,4)], file = paste('trainingFiles.impc/fcs/', sample, '.live.parent.fcs.rds', sep = ''))

        # Gating Live to obtain the Lymphocyte population. Plotting FSC-A_SSC-A
        lymph.flowD.temp <- flowDensity(live, channels = c(1,4), position = c(T, F), gates = c(gthres[5], gthres[6]))
        lymph.flowD <- flowDensity(lymph.flowD.temp, channels = c(1,4), position = c(F,F), gates = c(gthres[4], gthres[6]))
        lymph.flowD@proportion <- (lymph.flowD@cell.count/live.flowD@cell.count)*100
        lymph <- getflowFrame(lymph.flowD)

        parentFrames[["Lymphocytes"]] <- live
        sampleTrainGates[["Lymphocytes"]] <- TrainingGate(
                     parentName = "Live",
                     channelA = 1,
                     channelB = 4,
                     thresholdALow = gthres[5],
                     thresholdAHigh = gthres[4],
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[6],
                     negate = F,
                     gateAssignments = 1:nrow(live) %in% lymph.flowD@index
                     )

        saveRDS(live@exprs[,c(1,4)], file = paste('trainingFiles.impc/fcs/', sample, '.lymphocytes.parent.fcs.rds', sep = ''))

        # Gating Lymphocytes to obtain the CD45+ population. Plotting CD45_CD43 (CD45+ population)
        lymph.temp <-rotateData(lymph,c(14,8),theta = pi/6.50)$data

        cd45.flowD.temp <- flowDensity(lymph.temp, channels = c(14,8), position = c(T,F), gates = c(gthres[7], gthres[9]))
        cd45.flowD <- flowDensity(cd45.flowD.temp, channels = c(14,8), position = c(F,F), gates = c(gthres[8], gthres[9]))

        cd45.flowD@filter <-  rotateData(cd45.flowD@filter,c(14,8),theta = -pi/6.50)$data
        cd45.flowD@flow.frame <- rotateData(getflowFrame(cd45.flowD),c(14,8),theta = -pi/6.50)$data
        cd45.flowD@proportion <- (cd45.flowD@cell.count/lymph.flowD@cell.count)*100
        cd45 <- getflowFrame(cd45.flowD)

        parentFrames[["CD45"]] <- lymph.temp
        sampleTrainGates[["CD45"]] <- TrainingGate(
                     parentName = "Lymphocytes",
                     channelA = 14,
                     channelB = 8,
                     thresholdALow = gthres[7],
                     thresholdAHigh = gthres[8],
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[9],
                     negate = F,
                     gateAssignments = 1:nrow(lymph) %in% cd45.flowD@index
                     )

        saveRDS(lymph.temp@exprs[,c(14,8)], file = paste('trainingFiles.impc/fcs/', sample, '.cd45.parent.fcs.rds', sep = ''))

        # Gating CD45+ to obtain Granulocyte Pre and NOT(Granulocyte Pre) populations. Plotting GR1_CD43
        NOT.granulocyte.flowD <- flowDensity(cd45, channels = c(10, 8), position = c(F,NA), gates = c(gthres[10], gthres[9]))
        NOT.granulocyte <- getflowFrame(NOT.granulocyte.flowD)

        parentFrames[["NOT(Granulocyte Pre)"]] <- cd45
        sampleTrainGates[["NOT(Granulocyte Pre)"]] <- TrainingGate(
                     parentName = "CD45",
                     channelA = 10,
                     channelB = 8,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[10],
                     thresholdBLow = NaN,
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(cd45) %in% NOT.granulocyte.flowD@index
                     )

        granulocyte <- cd45
        granulocyte@exprs <- granulocyte@exprs[-NOT.granulocyte.flowD@index,]
        
        saveRDS(cd45@exprs[,c(10,8)], file = paste('trainingFiles.impc/fcs/', sample, '.notgranulocytepre.parent.fcs.rds', sep = ''))

        parentFrames[["Granulocyte Pre"]] <- cd45
        sampleTrainGates[["Granulocyte Pre"]] <- TrainingGate(
                     parentName = "CD45",
                     channelA = 10,
                     channelB = 8,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[10],
                     thresholdBLow = NaN,
                     thresholdBHigh = NaN,
                     negate = T,
                     gateAssignments = 1:nrow(cd45) %in% setdiff(1:nrow(cd45), NOT.granulocyte.flowD@index)
                     )

        saveRDS(cd45@exprs[,c(10,8)], file = paste('trainingFiles.impc/fcs/', sample, '.granulocytepre.parent.fcs.rds', sep = ''))

        # Gating NOT(Granulocyte Pre) population to obtain CD3 Tcells and NOT(CD3 Tcells). Plotting B220_CD3
        cd3.Tcell.flowD <- flowDensity(NOT.granulocyte, channels = c(18,16), position = c(F, T), gates = c(gthres[11], gthres[12]))

        parentFrames[["CD3 Tcell"]] <- NOT.granulocyte
        sampleTrainGates[["CD3 Tcell"]] <- TrainingGate(
                     parentName = "NOT(Granulocyte Pre)",
                     channelA = 18,
                     channelB = 16,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[11],
                     thresholdBLow = gthres[12],
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(NOT.granulocyte) %in% cd3.Tcell.flowD@index
                     )

        NOT.cd3.Tcell <- NOT.granulocyte
        NOT.cd3.Tcell@exprs <- NOT.cd3.Tcell@exprs[-cd3.Tcell.flowD@index,]        

        saveRDS(NOT.granulocyte@exprs[,c(18,16)], file = paste('trainingFiles.impc/fcs/', sample, '.cd3tcell.parent.fcs.rds', sep = ''))

        parentFrames[["NOT(CD3 Tcell)"]] <- NOT.granulocyte
        sampleTrainGates[["NOT(CD3 Tcell)"]] <- TrainingGate(
                     parentName = "NOT(Granulocyte Pre)",
                     channelA = 18,
                     channelB = 16,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[11],
                     thresholdBLow = gthres[12],
                     thresholdBHigh = NaN,
                     negate = T,
                     gateAssignments = 1:nrow(NOT.granulocyte) %in% setdiff(1:nrow(NOT.granulocyte), cd3.Tcell.flowD@index)
                     )        

        saveRDS(NOT.granulocyte@exprs[,c(18,16)], file = paste('trainingFiles.impc/fcs/', sample, '.notcd3tcell.parent.fcs.rds', sep = ''))

        # Gating NOT(CD3 Tcell) population to obtain Plasma and NOT Plasma. Plotting CD138_B220
        plasma.flowD <- flowDensity(NOT.cd3.Tcell, channels = c(15,18), position = c(T, F), gates = c(gthres[13], gthres[11]))        

        parentFrames[["Plasma"]] <- NOT.cd3.Tcell
        sampleTrainGates[["Plasma"]] <- TrainingGate(
                     parentName = "NOT(CD3 Tcell)",
                     channelA = 15,
                     channelB = 18,
                     thresholdALow = gthres[13],
                     thresholdAHigh = NaN,
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[11],
                     negate = F,
                     gateAssignments = 1:nrow(NOT.cd3.Tcell) %in% plasma.flowD@index
                     )

        NOT.plasma <- NOT.cd3.Tcell
        NOT.plasma@exprs <-  NOT.plasma@exprs[-plasma.flowD@index,]

        saveRDS(NOT.cd3.Tcell@exprs[,c(15,18)], file = paste('trainingFiles.impc/fcs/', sample, '.plasma.parent.fcs.rds', sep = ''))

        parentFrames[["NOT(Plasma)"]] <- NOT.cd3.Tcell
        sampleTrainGates[["NOT(Plasma)"]] <- TrainingGate(
                     parentName = "NOT(CD3 Tcell)",
                     channelA = 15,
                     channelB = 18,
                     thresholdALow = gthres[13],
                     thresholdAHigh = NaN,
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[11],
                     negate = T,
                     gateAssignments = 1:nrow(NOT.cd3.Tcell) %in% setdiff(1:nrow(NOT.cd3.Tcell), plasma.flowD@index)
                     )

        saveRDS(NOT.cd3.Tcell@exprs[,c(15,18)], file = paste('trainingFiles.impc/fcs/', sample, '.notplasma.parent.fcs.rds', sep = ''))

        # Gating NOT Plasma population to obtain Myeloid Pre and Bcells. Plotting B220_CD11b
        Bcell.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(T, NA), gates = c(gthres[11], gthres[14]))
        Bcell <- getflowFrame(Bcell.flowD)

        parentFrames[["Bcell"]] <- NOT.plasma
        sampleTrainGates[["Bcell"]] <- TrainingGate(
                     parentName = "NOT(Plasma)",
                     channelA = 18,
                     channelB = 13,
                     thresholdALow = gthres[11],
                     thresholdAHigh = NaN,
                     thresholdBLow = NaN,
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(NOT.plasma) %in% Bcell.flowD@index
                     )

        saveRDS(NOT.plasma@exprs[,c(18,13)], file = paste('trainingFiles.impc/fcs/', sample, '.bcell.parent.fcs.rds', sep = ''))

        myeloid.flowD <- flowDensity(NOT.plasma, channels = c(18,13), position = c(F, NA), gates = c(gthres[11], gthres[14]))

        parentFrames[["Myeloid"]] <- NOT.plasma
        sampleTrainGates[["Myeloid"]] <- TrainingGate(
                     parentName = "NOT(Plasma)",
                     channelA = 18,
                     channelB = 13,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[11],
                     thresholdBLow = NaN,
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(NOT.plasma) %in% myeloid.flowD@index
                     )

        saveRDS(NOT.plasma@exprs[,c(18,13)], file = paste('trainingFiles.impc/fcs/', sample, '.myeloid.parent.fcs.rds', sep = ''))

        # Gating Bcells population to obtain CD43+ and CD43-. Plotting B220_CD43
        cd43.pos.flowD <- flowDensity(Bcell, channels = c(18,8), position = c(NA, T), gates = c(gthres[11], gthres[15]))
        cd43.pos <- getflowFrame(cd43.pos.flowD)

        parentFrames[["CD43+"]] <- Bcell
        sampleTrainGates[["CD43+"]] <- TrainingGate(
                     parentName = "Bcell",
                     channelA = 18,
                     channelB = 8,
                     thresholdALow = NaN,
                     thresholdAHigh = NaN,
                     thresholdBLow = gthres[15],
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(Bcell) %in% cd43.pos.flowD@index
                     )

        saveRDS(Bcell@exprs[,c(18,8)], file = paste('trainingFiles.impc/fcs/', sample, '.cd43p.parent.fcs.rds', sep = ''))

        cd43.neg <- Bcell
        cd43.neg@exprs <- cd43.neg@exprs[-cd43.pos.flowD@index,]

        parentFrames[["CD43-"]] <- Bcell
        sampleTrainGates[["CD43-"]] <- TrainingGate(
                     parentName = "Bcell",
                     channelA = 18,
                     channelB = 8,
                     thresholdALow = NaN,
                     thresholdAHigh = NaN,
                     thresholdBLow = gthres[15],
                     thresholdBHigh = NaN,
                     negate = T,
                     gateAssignments = 1:nrow(Bcell) %in% setdiff(1:nrow(Bcell), cd43.pos.flowD@index)
                     )

        saveRDS(Bcell@exprs[,c(18,8)], file = paste('trainingFiles.impc/fcs/', sample, '.cd43n.parent.fcs.rds', sep = ''))

        # Gating CD43+ population to obtain HFA, HFB, and HFC. Plotting CD24_BP1
        HFA.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(F,F), gates = c(gthres[16], gthres[19]))

        parentFrames[["HFA"]] <- cd43.pos
        sampleTrainGates[["HFA"]] <- TrainingGate(
                     parentName = "CD43+",
                     channelA = 9,
                     channelB = 17,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[16],
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[19],
                     negate = F,
                     gateAssignments = 1:nrow(cd43.pos) %in% HFA.flowD@index
                     )

        saveRDS(cd43.pos@exprs[,c(9,17)], file = paste('trainingFiles.impc/fcs/', sample, '.hfa.parent.fcs.rds', sep = ''))

        HFB.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(T,F), gates = c(gthres[16], gthres[19]))

        parentFrames[["HFB"]] <- cd43.pos
        sampleTrainGates[["HFB"]] <- TrainingGate(
                     parentName = "CD43+",
                     channelA = 9,
                     channelB = 17,
                     thresholdALow = gthres[16],
                     thresholdAHigh = NaN,
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[19],
                     negate = F,
                     gateAssignments = 1:nrow(cd43.pos) %in% HFB.flowD@index
                     )

        saveRDS(cd43.pos@exprs[,c(9,17)], file = paste('trainingFiles.impc/fcs/', sample, '.hfb.parent.fcs.rds', sep = ''))

        HFC.flowD <- flowDensity(cd43.pos, channels = c(9, 17), position = c(T,T), gates = c(gthres[16], gthres[19]))

        parentFrames[["HFC"]] <- cd43.pos
        sampleTrainGates[["HFC"]] <- TrainingGate(
                     parentName = "CD43+",
                     channelA = 9,
                     channelB = 17,
                     thresholdALow = gthres[16],
                     thresholdAHigh = NaN,
                     thresholdBLow = gthres[19],
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(cd43.pos) %in% HFC.flowD@index
                     )

        saveRDS(cd43.pos@exprs[,c(9,17)], file = paste('trainingFiles.impc/fcs/', sample, '.hfc.parent.fcs.rds', sep = ''))

        # Gating CD43- population to obtain HFD, HFE, and HFF. Plotting IgM_IgD
        HFD.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(F,F), gates = c(gthres[17], gthres[18]))

        parentFrames[["HFD"]] <- cd43.neg
        sampleTrainGates[["HFD"]] <- TrainingGate(
                     parentName = "CD43-",
                     channelA = 12,
                     channelB = 7,
                     thresholdALow = NaN,
                     thresholdAHigh = gthres[17],
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[18],
                     negate = F,
                     gateAssignments = 1:nrow(cd43.neg) %in% HFD.flowD@index
                     )

        saveRDS(cd43.neg@exprs[,c(12,7)], file = paste('trainingFiles.impc/fcs/', sample, '.hfd.parent.fcs.rds', sep = ''))

        HFE.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(T,F), gates = c(gthres[17], gthres[18]))

        parentFrames[["HFE"]] <- cd43.neg
        sampleTrainGates[["HFE"]] <- TrainingGate(
                     parentName = "CD43-",
                     channelA = 12,
                     channelB = 7,
                     thresholdALow = gthres[17],
                     thresholdAHigh = NaN,
                     thresholdBLow = NaN,
                     thresholdBHigh = gthres[18],
                     negate = F,
                     gateAssignments = 1:nrow(cd43.neg) %in% HFE.flowD@index
                     )

        saveRDS(cd43.neg@exprs[,c(12,7)], file = paste('trainingFiles.impc/fcs/', sample, '.hfe.parent.fcs.rds', sep = ''))

        HFF.flowD <- flowDensity(cd43.neg, channels = c(12,7), position = c(NA,T), gates = c(gthres[17], gthres[18]))

        parentFrames[["HFF"]] <- cd43.neg
        sampleTrainGates[["HFF"]] <- TrainingGate(
                     parentName = "CD43-",
                     channelA = 12,
                     channelB = 7,
                     thresholdALow = NaN,
                     thresholdAHigh = NaN,
                     thresholdBLow = gthres[18],
                     thresholdBHigh = NaN,
                     negate = F,
                     gateAssignments = 1:nrow(cd43.neg) %in% HFF.flowD@index
                     )

        saveRDS(cd43.neg@exprs[,c(12,7)], file = paste('trainingFiles.impc/fcs/', sample, '.hff.parent.fcs.rds', sep = ''))

        for (pop in names(sampleTrainGates))
        {
            tg <- sampleTrainGates[[pop]]
            pf <- parentFrames[[pop]]

            sampleTrainGates[[pop]]@densitiesA[[as.character(numFeatures)]] <- estimateDensity(pf@exprs[, tg@channelA], numFeatures)
            sampleTrainGates[[pop]]@densitiesB[[as.character(numFeatures)]] <- estimateDensity(pf@exprs[, tg@channelB], numFeatures)
        }

        saveRDS(sampleTrainGates, file = paste('trainingFiles.impc/', sample, '.rds', sep = ''))

    },
    error = function(err)
    {
        print(err)
    })

})

tr <- list()

tr[['singlets']] <- readTrainFiles('trainingFiles.impc', 'Singlets', numFeatures)
tr[['live']] <- readTrainFiles('trainingFiles.impc', 'Live', numFeatures)
tr[['lymphocytes']] <- readTrainFiles('trainingFiles.impc', 'Lymphocytes', numFeatures)
tr[['cd45']] <- readTrainFiles('trainingFiles.impc', 'CD45', numFeatures)
tr[['notgranulocytepre']] <- readTrainFiles('trainingFiles.impc', 'NOT(Granulocyte Pre)', numFeatures)
tr[['granulocytepre']] <- readTrainFiles('trainingFiles.impc', 'Granulocyte Pre', numFeatures)
tr[['cd3tcell']] <- readTrainFiles('trainingFiles.impc', 'CD3 Tcell', numFeatures)
tr[['notcd3tcell']] <- readTrainFiles('trainingFiles.impc', 'NOT(CD3 Tcell)', numFeatures)
tr[['plasma']] <- readTrainFiles('trainingFiles.impc', 'Plasma', numFeatures)
tr[['notplasma']] <- readTrainFiles('trainingFiles.impc', 'NOT(Plasma)', numFeatures)
tr[['bcell']] <- readTrainFiles('trainingFiles.impc', 'Bcell', numFeatures)
tr[['myeloid']] <- readTrainFiles('trainingFiles.impc', 'Myeloid', numFeatures)
tr[['cd43p']] <- readTrainFiles('trainingFiles.impc', 'CD43+', numFeatures)
tr[['cd43n']] <- readTrainFiles('trainingFiles.impc', 'CD43-', numFeatures)
tr[['hfa']] <- readTrainFiles('trainingFiles.impc', 'HFA', numFeatures)
tr[['hfb']] <- readTrainFiles('trainingFiles.impc', 'HFB', numFeatures)
tr[['hfc']] <- readTrainFiles('trainingFiles.impc', 'HFC', numFeatures)
tr[['hfd']] <- readTrainFiles('trainingFiles.impc', 'HFD', numFeatures)
tr[['hfe']] <- readTrainFiles('trainingFiles.impc', 'HFE', numFeatures)
tr[['hff']] <- readTrainFiles('trainingFiles.impc', 'HFF', numFeatures)

save(list = c('tr'), file = 'trainingFiles.impc/tr.RData')

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

save(list = c('tr', 'dA', 'dB'), file = 'trainingFiles.impc/tr.dtwdists.RData')