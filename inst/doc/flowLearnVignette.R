## ------------------------------------------------------------------------
library(flowLearn)

## ------------------------------------------------------------------------
samples <- unique(flSampleDensdat@data$fcs)
populations <- unique(flSampleDensdat@data$population)
print(samples)
print(populations)

## ---- fig.width=7, fig.height=5------------------------------------------
# Select the first sample as an example.
sampleIdx <- 1
fcsName <- samples[[sampleIdx]]

# Use flFind() to filter the second channel density of the notplasma population with the given sample.
dd <- flFind(flSampleDensdat, paste0('population == "notplasma" & channelIdx == 2 & fcs == "', fcsName, '"'))

# Use flPlotDensThresh() to plot the selected density and threshold.
flPlotDensThresh(flGetDensity(dd), flGetGate(dd))

## ------------------------------------------------------------------------
# Extract all densities of the bcell population and the first channel.
dd <- flFind(flSampleDensdat, 'population == "bcell" & channelIdx == 1')

# Use flSelectPrototypes() to select one prototype.
protoIdx <- flSelectPrototypes(dd, 1)

print(paste0('Selected prototype index: ', protoIdx))

## ---- fig.width=7, fig.height=5------------------------------------------
# Use flPredictThresholds() to predict thresholds using a given prototype
# All other samples are aligned with the prototype and its threshold is transferred.
ddp <- flPredictThresholds(dd, protoIdx)
                       
# Plot predictions the prediction, exemplary for the third sample
# Using flAt(), the density at index 3 is extracted for display.
# flGetDensity() and flGetGate() extract its density and gate, respectively.
flPlotDensThresh(flGetDensity(flAt(ddp, 3)), flGetGate(flAt(dd, 3)), flGetGate(flAt(ddp, 3)))

## ---- fig.width=7, fig.height=5------------------------------------------
# Re-do the alignment for the purpose of visualization.
# The third density is aligned to the prototype density.
# Here 
d <- flDtwMain(flGetDensity(flAt(dd, 3)), flGetDensity(flAt(dd, protoIdx)))

# Use the dtw package's plotting routine.
plot(d, type = 'twoway', lwd = 2, offset = 0.001, match.indices = 200)

## ---- fig.width=7, fig.height=5------------------------------------------
# Calculate the F1 score for each sample
f1Scores <- sapply(samples, function(fcs) 
{

	# The function flEvalF1ScoreFCS() takes the predicted DensityData object, an FCS sample name, 
	# the predicted population name, true gate assignments and parent population expressions.
	# It gates the child population using the predicted gates and compares the result to 
	# the true population in terms of precision, recall and effectively F1 score.

	flEvalF1ScoreFCS(ddp, 
					 fcs, 
					 'bcell', 
					 flSampleBcellEvaluationData[[fcs]]$gateAssignments, 
					 flSampleBcellEvaluationData[[fcs]]$parentExprs, 
					 F)

})

print(f1Scores)

## ------------------------------------------------------------------------
densdat <- flInit(new('DensityData'))

## ------------------------------------------------------------------------
# Use flEstimateDensity() to obtain a kernel density estimate that is smooth.
# It is given a vector of cell measurements for one channel and the number of density features.
# In this vignette example, the provided flSampleFlowFrame contains only two channels. 

densA <- flEstimateDensity(flSampleFlowFrame@exprs[, 1], densdat@numFeatures)
densB <- flEstimateDensity(flSampleFlowFrame@exprs[, 1], densdat@numFeatures)

# Now add the densities to the `densdat` object.
densdat <- flAdd(densdat, flSampleFlowFrame@description$"$FIL", 'granulocytepre', 1, densA$x, densA$y)
densdat <- flAdd(densdat, flSampleFlowFrame@description$"$FIL", 'granulocytepre', 2, densB$x, densB$y)

## ------------------------------------------------------------------------
print(flSize(densdat) == 2)

