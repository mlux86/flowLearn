library(stringr)
library(flowLearn)
library(FlowSOM)
library(parallel)
library(plyr)
library(ggplot2)

source('helper_match_evaluate_multiple.R')

printf <- function(...) invisible(print(sprintf(...)))

bonemarrowPath <- '/mnt/data/BM_Panel/'
cd45Path <- paste0(bonemarrowPath, 'CleanFiles-CD45/')


load(paste0(bonemarrowPath, 'store.allFCS-Updated.Rdata'))
load(paste0(bonemarrowPath, 'all.gthres.Updated.Rdata'))

sampleMeta <- as.data.frame(store.allFCS, stringsAsFactors = F)
threshTable <- as.data.frame(all.gthres.Updated, stringsAsFactors = F)
gateMeta <- read.table('prepare_bonemarrow.tsv', sep = '\t', header = T, stringsAsFactors = F)

rownames(sampleMeta) <- sampleMeta$"FCS files"
rownames(threshTable) <- sampleMeta$"FCS files"

# generate gating infos

gatingInfos <- list()
for (i in 1:nrow(gateMeta))
{
	popName <- gateMeta[i, "population"]
	parentName <- gateMeta[i, "parent"]
	channelIndices <- c(gateMeta[i, "channelA"], gateMeta[i, "channelB"])
	gatingInfos[[popName]] <- new("GatingInfo", population = popName, parent = parentName, channels = channelIndices)	
}

# generate gating ground truth

filenames <- dir(cd45Path)

cl <- parallel::makeCluster(parallel::detectCores(), type = "FORK", outfile = "")
f1s <- parLapply(cl, filenames, function(fname) 
# f1s <- parLapply(cl, filenames[sample(length(filenames), 100)], function(fname) 
# f1s <- lapply(filenames, function(fname)
{
	barcode <- sub('.*(L[0-9]+).*', '\\1', fname)
	
	print(barcode)

	tryCatch(
	{

		meta <- subset(sampleMeta, Barcodes == barcode)
		thresh <- unlist(threshTable[grep(barcode, rownames(threshTable)), ])
		names(thresh) <- NULL

		load(paste0(cd45Path, fname))
        nCd45 <- nrow(cd45@exprs)

        # apply flowSOM on cd45 flow frame
        # start.time <- Sys.time()
        res.flowSOM <- FlowSOM(cd45, colsToUse = 7:19, nClus = 11)
        # end.time <- Sys.time()
        # time.taken <- end.time - start.time
        # print(paste0("FlowSOM time taken for one sample:", time.taken))

		frames <- list()
		frames$cd45 <- cd45	

        lstParentIndices <- list()
        lstParentNames <- list()

        popLabels <- character(nCd45)

		for (i in 1:nrow(gateMeta))
		{

			popName <- gateMeta[i, "population"]
			parentName <- gateMeta[i, "parent"]
			channelIndices <- c(gateMeta[i, "channelA"], gateMeta[i, "channelB"])
			thresholds <- list(thresholdALow = thresh[gateMeta[i, "gateALow"]], thresholdAHigh = thresh[gateMeta[i, "gateAHigh"]], thresholdBLow = thresh[gateMeta[i, "gateBLow"]], thresholdBHigh = thresh[gateMeta[i, "gateBHigh"]])

			parentFrame <- frames[[parentName]]

			parentIndices <- replicate(nrow(parentFrame@exprs), T)
			if (!is.na(thresholds$thresholdALow))
			{ 
				parentIndices <- parentIndices & (parentFrame@exprs[, channelIndices[1]] > thresholds$thresholdALow)
			}
			if (!is.na(thresholds$thresholdAHigh))
			{ 
				parentIndices <- parentIndices & (parentFrame@exprs[, channelIndices[1]] < thresholds$thresholdAHigh)
			}
			if (!is.na(thresholds$thresholdBLow))
			{ 
				parentIndices <- parentIndices & (parentFrame@exprs[, channelIndices[2]] > thresholds$thresholdBLow)
			}
			if (!is.na(thresholds$thresholdBHigh))
			{ 
				parentIndices <- parentIndices & (parentFrame@exprs[, channelIndices[2]] < thresholds$thresholdBHigh)
			}

			if (gateMeta[i, 'negate'])
			{
				parentIndices <- !parentIndices
			}

			popFrame <- parentFrame
			popFrame@exprs <- popFrame@exprs[parentIndices, ]
			frames[[popName]] <- popFrame

            # calculate root indices w.r.t. cd45 frame, not parent frame
            lstParentIndices[[popName]] <- parentIndices
            lstParentNames[[popName]] <- parentName

            cpLstParentIndices <- lstParentIndices
            popName_ <- popName
            while(lstParentNames[[popName_]] != 'cd45')
            {
                parentName <- lstParentNames[[popName_]]
                cpLstParentIndices[[parentName]][which(cpLstParentIndices[[parentName]])] <- cpLstParentIndices[[popName_]]
                popName_ <- parentName
            }

            rootIdx <- cpLstParentIndices[[popName_]]

            if(!(popName %in% c("notgranulocytepre", "notcd3tcell", "notplasma", "bcell", "cd43p", "cd43n")))
            {
                popLabels[rootIdx] <- popName
            }

		}

        lblCorrect <- as.factor(popLabels)
        popNamesCorrect <- levels(lblCorrect)
        levels(lblCorrect) <- 1:length(popNamesCorrect)
        lblCorrect <- as.numeric(lblCorrect)

        lblFlowSOM <- as.numeric(res.flowSOM[[2]][res.flowSOM[[1]]$map$mapping[,1]])

        z <- helper_match_evaluate_multiple(lblFlowSOM, lblCorrect)

        f1 <- z$F1

        popNamesCorrect[popNamesCorrect == ""] <- "ungated"
        names(f1) <- popNamesCorrect

        return(f1)

	}, error = function(e)
	{
		print(paste0('Error in file "', fname, '":'))
		print(e)
		return(NA)
	})

})

# l <- laply(f1s, length)
# f1s <- as.data.frame(do.call('rbind', f1s[l == 11]))

# colnames(f1s) <- c("ungated", "CD3 T-cell", "Granulocyte Pre", "HFA", "HFB", "HFC", "HFD", "HFE", "HFF", "Myeloid", "Plasma")

# print(paste0("Samples with wrong number of clusters: ", sum(l != 11)))

# save(f1s, file = 'eval_flowsom_f1.RData')

# f1s <- f1s[-1]

# p <- ggplot(stack(f1s), aes(x = ind, y = values)) +
#     stat_boxplot(geom = 'errorbar', width = 0.25) +
#     geom_boxplot(width = 0.3, outlier.shape = 20, outlier.size = 0.1) +
#     scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.05)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     xlab('Population') +
#     ylab('F1-score')

# print(p)

# ggsave('results/eval_flowsom_f1_impc.bonemarrow.png')    
# ggsave('results/eval_flowsom_f1_impc.bonemarrow.eps')    
