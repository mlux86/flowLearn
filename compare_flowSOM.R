library(stringr)
library(flowLearn)
library(FlowSOM)
library(parallel)

source('helper_match_evaluate_multiple.R')

printf <- function(...) invisible(print(sprintf(...)))

bonemarrowPath <- '/mnt/data/BM_Panel/'
cd45Path <- paste0(bonemarrowPath, 'CleanFiles-CD45/')
p <- '/home/mlux/flowlearn.traindat.impc.bonemarrow/'


dir.create(p, showWarnings = FALSE)
dir.create(paste0(p, 'eval'), showWarnings = FALSE)

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
        res.flowSOM <- FlowSOM(cd45, colsToUse = 7:19, nClus = 11)

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
        # names(f1) <- popNamesCorrect
        # print(f1)
        # print(paste0(z$mean_F1, "  ----------"))

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

save(f1s, file = 'f1s.RData')

