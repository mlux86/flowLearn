library(stringr)
library(parallel)
library(flowLearn)
library(flowCore)

printf <- function(...) invisible(print(sprintf(...)))

micePath <- '/mnt/data/BM_Panel/'
cd45Path <- paste0(micePath, 'CleanFiles-CD45/')
p <- '/home/mlux/flowlearn.traindat.mice/'


dir.create(p, showWarnings = FALSE)
dir.create(paste0(p, 'eval'), showWarnings = FALSE)

load(paste0(micePath, 'store.allFCS-Updated.Rdata'))
load(paste0(micePath, 'all.gthres.Updated.Rdata'))

sampleMeta <- as.data.frame(store.allFCS, stringsAsFactors = F)
threshTable <- as.data.frame(all.gthres.Updated, stringsAsFactors = F)
gateMeta <- read.table('prepare_mice.tsv', sep = '\t', header = T, stringsAsFactors = F)

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
densdat.lst <- parLapply(cl, filenames, function(fname)
{
	barcode <- sub('.*(L[0-9]+).*', '\\1', fname)

	print(barcode)

	tryCatch(
	{

		meta <- subset(sampleMeta, Barcodes == barcode)
		thresh <- unlist(threshTable[grep(barcode, rownames(threshTable)), ])
		names(thresh) <- NULL

		# load(paste0(cd45Path, fname))
		cd4 <- read.FCS(cd45Path)

		frames <- list()
		frames$cd45 <- cd45

		densdat <- new('DensityData')

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

			densA <- flEstimateDensity(parentFrame@exprs[, channelIndices[1]], 512)
			densB <- flEstimateDensity(parentFrame@exprs[, channelIndices[2]], 512)

			densdat <- flAdd(densdat, meta$"FCS files", popName, 1, densA$x, densA$y, thresholds$thresholdALow, thresholds$thresholdAHigh)
			densdat <- flAdd(densdat, meta$"FCS files", popName, 2, densB$x, densB$y, thresholds$thresholdBLow, thresholds$thresholdBHigh)

			evaldat <- list(parentFrame = parentFrame@exprs[, channelIndices], indices = parentIndices, negate = gateMeta[i, 'negate'])

			saveRDS(evaldat, file = paste0(p, 'eval/', meta$"FCS files", '.', popName, '.eval.rds'))

		}

		return(densdat)

	}, error = function(e)
	{
		print(paste0('Error in file "', fname, '":'))
		print(e)
		return(NA)
	})

})

parallel::stopCluster(cl)

densdat <- Reduce(function(x, y)
{
	if(is.na(y))
	{
		return(x)
	}
	return(flConcat(x, y))
}, densdat.lst, flInit(new('DensityData')))
save(densdat, gatingInfos, sampleMeta, file = paste0(p, 'train_data.RData'))

