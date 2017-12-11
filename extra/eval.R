library(ggplot2)
library(gridExtra)
library(Rtsne)
library(plyr)
library(flowLearn)

printf <- function(...) invisible(cat(sprintf(...)))

flEvalF1ScorePopulation <- function(densdat, population)
{
    fcsFiles <- as.character(unique(flData(densdat)$fcs))
    df <- ldply(fcsFiles, function(fcs)
                {
                    return(flEvalF1ScoreFCS(densdat, fcs, population))
                })
    rownames(df) <- fcsFiles
    return(df)
}

flEvalDataset <- function(datasetName, numProtoPerChannel, traindatFolderPrefix = '~/flowlearn.traindat/flowlearn.traindat.')
{
    flowLearnEvaluationFolder <<- paste0(traindatFolderPrefix, datasetName, '/')

    load(paste0(flowLearnEvaluationFolder, 'train_data.RData'))

    if(file.exists('popmap.csv'))
    {
        pm <- read.csv('popmap.csv', stringsAsFactors = F)
        rownames(pm) <- pm$normalized
    }

    fcsFiles <- unique(flData(densdat)$fcs)
    nSamples <- length(fcsFiles)
    populations <- as.character(unique(flData(densdat)$population))

    f1s <- matrix(NaN, nSamples, length(populations)) # F1 scores
    tp <- matrix(NaN, nSamples, length(populations)) # true proportions
    pp <- matrix(NaN, nSamples, length(populations)) # predicted proportions

    rownames(f1s) <- fcsFiles
    rownames(tp) <- fcsFiles
    rownames(pp) <- fcsFiles

    for (i in 1:length(populations))
    {

        population <- populations[i]


        tryCatch({

            # start.time <- Sys.time()

            # predict
            dd1 <- flFind(densdat, population = population, channelIdx = 1)
            dd2 <- flFind(densdat, population = population, channelIdx = 2)
            protoIdx1 <- flSelectPrototypes(dd1, numProtoPerChannel)
            protoIdx2 <- flSelectPrototypes(dd2, numProtoPerChannel)
            ddp1 <- flPredictThresholds(dd1, protoIdx1)
            ddp2 <- flPredictThresholds(dd2, protoIdx2)

            # end.time <- Sys.time()
            # time.taken <- end.time - start.time
            # print(paste0("FlowLearn time taken for one population: ", time.taken))

            # evaluate
            protoIdx <- union(protoIdx1, protoIdx2)
            testIdx <- setdiff(1:nSamples, protoIdx)
            ddp1test <- flAt(ddp1, testIdx)
            ddp2test <- flAt(ddp2, testIdx)
            ddp <- flConcat(ddp1test, ddp2test)

            e <- flEvalF1ScorePopulation(ddp, population)
            f1s[testIdx, i] <- e$f1
            tp[testIdx, i] <- e$trueProportion
            pp[testIdx, i] <- e$predictedProportion

        }, error = function(e)
        {
            print(paste0("Error in population ", population))
            print(e)
        })

    }

    for(ip in 1:length(populations))
    {
        normpop <- populations[[ip]]
        populations[[ip]] <- pm[normpop,]$original
    }

    dfEval <- as.data.frame(f1s)
    colnames(dfEval) <- populations
    dfEvalTp <- as.data.frame(tp)
    colnames(dfEvalTp) <- populations
    dfEvalPp <- as.data.frame(pp)
    colnames(dfEvalPp) <- populations

    nanidx <- which(sapply(dfEval, function(x) { sum(is.nan(x)) == nrow(dfEval) }))
    if(length(nanidx) > 0)
    {
        dfEval <- dfEval[, -nanidx]
        dfEvalTp <- dfEvalTp[, -nanidx]
        dfEvalPp <- dfEvalPp[, -nanidx]
    }

    save(dfEval, dfEvalTp, dfEvalPp, numProtoPerChannel, file = sprintf('results/eval_%s_%02d.RData', datasetName, numProtoPerChannel))

    p <- ggplot(stack(dfEval), aes(x = ind, y = values)) +
        stat_boxplot(geom = 'errorbar', width = 0.25) +
        geom_boxplot(width = 0.3, outlier.shape = 20, outlier.size = 0.1) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.05)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab('Population') +
        ylab('F1-score')

    print(p)

    ggsave(sprintf('results/eval_f1_%s_%02d.png', datasetName, numProtoPerChannel))
    ggsave(sprintf('results/eval_f1_%s_%02d.eps', datasetName, numProtoPerChannel))

    s1 <- stack(dfEvalTp)
    s1$proportionType = 'true'
    s2 <- stack(dfEvalPp)
    s2$proportionType = 'predicted'
    s <- rbind(s1, s2)

    p2 <- ggplot(s, aes(x = ind, y = values, fill = proportionType)) +
        stat_boxplot(geom = 'errorbar', width = 0.25) +
        geom_boxplot(width = 0.5, outlier.shape = 20, outlier.size = 0.1, notch = T) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.05)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab('Population') +
        ylab('Proportion')

    print(p2)

    ggsave(sprintf('results/eval_proportion_%s_%02d.png', datasetName, numProtoPerChannel))
    ggsave(sprintf('results/eval_proportion_%s_%02d.eps', datasetName, numProtoPerChannel))
}

flPlotPredictions <- function(densdat, pop, chan, numProto, oneplot = T, save = F)
{
    nSamples <- length(unique(flData(densdat)$fcs))

    dd <- flFind(densdat, population = pop, channelIdx = chan)
    protoIdx <- flSelectPrototypes(dd, numProto)
    ddp <- flPredictThresholds(dd, protoIdx)
    D <- as.matrix(dist(flGetDensity(dd)$y, method = "manhattan"))

    if(!oneplot)
    {
        graphics.off()
        dev.new()
        dev.new()
    }

    for (i in 1:nSamples)
    {
        if(oneplot)
        {
            par(mfrow=c(2,1))
        } else {
            dev.set(dev.list()[[1]])
        }
        flPlotDensThresh(flGetDensity(flAt(dd, i)), flGetGate(flAt(dd, i)), flGetGate(flAt(ddp, i)))
        if(save)
        {
            png('/tmp/prediction.png', width = 2000, height = 1500, res = 300)
            flPlotDensThresh(flGetDensity(flAt(dd, i)), flGetGate(flAt(dd, i)), flGetGate(flAt(ddp, i)))
            dev.off()
        }
        j <- order(D[protoIdx, i])[1]
        d <- flDtwMain(flGetDensity(flAt(dd, i)), flGetDensity(flAt(dd, protoIdx[j])))
        if(!oneplot)
        {
            dev.set(dev.list()[[2]])
        }
        plot(d, type = 'twoway', lwd = 2, offset = 0.001, match.indices = 100)
        if(save)
        {
            png('/tmp/alignment.png', width = 2000, height = 1500, res = 300)
            plot(d, type = 'twoway', lwd = 2, offset = 0.001, match.indices = 100)
            dev.off()
        }
        readline()
    }
}

flPlotPopulation <- function(densdat, pop)
{
    nSamples <- length(unique(flData(densdat)$fcs))

    dd1 <- flFind(densdat, population = pop, channelIdx = 1)
    dd2 <- flFind(densdat, population = pop, channelIdx = 2)

    ysne1 <- Rtsne(flGetDensity(dd1)$y, perplexity = log(nSamples^2))$Y
    ysne2 <- Rtsne(flGetDensity(dd2)$y, perplexity = log(nSamples^2))$Y
    ypca1 <- prcomp(flGetDensity(dd1)$y)$x[,1:2]
    ypca2 <- prcomp(flGetDensity(dd2)$y)$x[,1:2]

    df1 <- data.frame(sne1 = ysne1[,1], sne2 = ysne1[,2], pca1 = ypca1[,1], pca2 = ypca1[,2])
    df2 <- data.frame(sne1 = ysne2[,1], sne2 = ysne2[,2], pca1 = ypca2[,1], pca2 = ypca2[,2])

    p1 <- ggplot(df1, aes(x=pca1, y=pca2)) + geom_point(size = 3)
    p2 <- ggplot(df1, aes(x=sne1, y=sne2)) + geom_point(size = 3)
    p3 <- ggplot(df2, aes(x=pca1, y=pca2)) + geom_point(size = 3)
    p4 <- ggplot(df2, aes(x=sne1, y=sne2)) + geom_point(size = 3)
    plist <- list(p1, p2, p3, p4)
    do.call("grid.arrange", c(plist, ncol = 2))
}
