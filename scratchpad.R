source('flowLearn.R')

load('trainingFiles/tr.dtwdists.RData')

numTrain <- 2

population <- 'CD64++CD16+'
# population <- 'CD14+CD16-'

popId <- normalizePopulationName(population)

tr <- tr[[popId]]
dA <- dA[[popId]]
dB <- dB[[popId]]

n <- nrow(tr@densYchanA)

protoIdx <- 4

# selectedPrototypes <- selectPrototypes(tr, dA, dB, numTrain)
selectedPrototypes <- selectFixedPrototypes(tr, dA, dB, c(protoIdx), c(protoIdx))

predictedThresholds <- predictThresholds(tr, selectedPrototypes)

evaluatePerformance(tr, population, selectedPrototypes$testIdx, predictedThresholds$threshA, predictedThresholds$threshB)





# for (z in 1:length(selectedPrototypes$testIdx))
# {

# 	i <- selectedPrototypes$testIdx[z]

# 	protoIdx <- selectedPrototypes$trainIdxA[selectedPrototypes$labelsA[i]]

# 	print('train')
# 	print(tr@samples[[protoIdx]])
# 	print('test')
# 	print(tr@samples[[i]])

# 	dtwObj <- myDtw(tr@densYchanA[i,], tr@densYchanA[protoIdx,], k = T)

# 	plot(dtwObj, type = 'twoway', offset = 3, match.indices = 100)

# 	dev.new()
# 	plotDensThresh(tr@densXchanA[protoIdx,], tr@densYchanA[protoIdx,], tr@threshA[protoIdx,])

# 	dev.new()
# 	plotDensThresh(tr@densXchanA[i,], tr@densYchanA[i,], tr@threshA[i,], predictedThresholds$threshA[z,])

# 	readline()

# }


for (z in 1:length(selectedPrototypes$testIdx))
{

	i <- selectedPrototypes$testIdx[z]

	protoIdx <- selectedPrototypes$trainIdxB[selectedPrototypes$labelsB[i]]

	print('train')
	print(tr@samples[[protoIdx]])
	print('test')
	print(tr@samples[[i]])

	dtwObj <- myDtw(tr@densYchanB[i,], tr@densYchanB[protoIdx,], k = T);	plot(dtwObj, type = 'twoway', offset = 3, match.indices = 100)

	# dev.new()
	# plotDensThresh(tr@densXchanB[protoIdx,], tr@densYchanB[protoIdx,], tr@threshB[protoIdx,])

	# dev.new()
	# plotDensThresh(tr@densXchanB[i,], tr@densYchanB[i,], tr@threshB[i,], predictedThresholds$threshB[z,])

	# readline()

}


