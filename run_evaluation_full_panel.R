library(flowLearn)

source('evaluation.R')

preloaded <- new.env()
load('trainingFiles/tr.dtwdists.RData', preloaded)
n <- length(preloaded$tr$cd45@samples) # we could've chosen any population here. It's just about the number of samples

# populations for the OneStudy data set

# populations <- list(
#                     'CD45',
#                     'PBMC',
#                     'Tcell',
#                     'Bcell',
#                     'DN',
#                     'CD4',
#                     'CD8',
#                     'CD14',
#                     'CD14++CD16-',
#                     'CD14++CD16+',
#                     'CD14+CD16-',
#                     'CD14+CD16+',
#                     'NK',
#                     'NKT',
#                     'CD56Bright',
#                     'CD64++CD16+',
#                     'CD56BrightCD16-'
#                     )


# populations for the IMPC bone marrow data set

populations <- list(
                    'Singlets',
                    'Live',
                    'Lymphocytes',
                    'CD45',
                    'NOT(Granulocyte Pre)',
                    'Granulocyte Pre',
                    'CD3 Tcell',
                    'NOT(CD3 Tcell)',
                    'Plasma',
                    'NOT(Plasma)',
                    'Bcell',
                    'Myeloid',
                    'CD43+',
                    'CD43-',
                    'HFA',
                    'HFB',
                    'HFC',
                    'HFD',
                    'HFE',
                    'HFF'
                    )


# how many prototypes (i.e. 1 per channel means numTrain <- 2)
numTrain <- 2

# how many runs? i.e how often to chose random prototypes?
sampleSize <- 10

meanF1 <- matrix(NaN, sampleSize, length(populations))
medianF1 <- matrix(NaN, sampleSize, length(populations))

for (i in 1:length(populations))
{

    population <- populations[[i]]

    for (j in 1:sampleSize)
    {

        perf <- evaluatePopulation(population, numTrain, seed = 42 + j, preloaded = preloaded)
        meanF1[j, i] <- perf$meanPerf[3]
        medianF1[j, i] <- perf$medianPerf[3]

    }

}

colnames(meanF1) <- populations
dfMeanF1 <- as.data.frame(meanF1)

colnames(medianF1) <- populations
dfMedianF1 <- as.data.frame(medianF1)

# proportions <- unlist(lapply(populations, function(pop) getMeanProportion(preloaded$tr, pop)))
# if (exists('proportions'))
# {
#     populationsAndProportions <- lapply(1:length(populations), function(i) sprintf('%s (%.2f%%)', populations[[i]], proportions[[i]]*100))
#     colnames(dfMedianF1) <- populationsAndProportions
#     colnames(dfMeanF1) <- populationsAndProportions
# }

p <- ggplot(stack(dfMedianF1), aes(x = ind, y = values)) + 
    geom_boxplot() + 
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.05)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste0('median f1 scores over ', sampleSize, ' different runs, each with\n', numTrain / 2, ' random prototype(s) per channel, tested on ', n - numTrain, ' samples')) + 
    xlab('Population and proportion') + 
    ylab('median(f1-score)')

print(p)