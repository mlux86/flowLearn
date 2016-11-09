library(flowLearn)

library(ggplot2)
library(Rtsne)
library(gridExtra)
library(stringr)

preloaded <- new.env()
load('trainingFiles.impc.new/tr.dtwdists.RData', preloaded)

population <- 'Singlets'
popId <- normalizePopulationName(population)

tr <- preloaded$tr[[popId]]
dA <- preloaded$dA[[popId]]
dB <- preloaded$dB[[popId]]

n <- nrow(tr@densYchanA)

# m <- as.dist(dA); 
m <- dB; 

pca <- princomp(m)
tsne <- Rtsne(m, dims = 2, perplexity = round(log(nrow(m))^2)); 
df <- data.frame(sne1 = tsne$Y[, 1], sne2 = tsne$Y[, 2], pca1 = pca$scores[, 1], pca2 = pca$scores[, 2], labels <- 1:n)

pPca <- ggplot(df, aes(x = pca1, y = pca2, colour = labels)) + geom_point(size = 2)
pSne <- ggplot(df, aes(x = sne1, y = sne2, colour = labels)) + geom_point(size = 2)

pdf(paste0("/tmp/pdfs/", popId, "_channelB.pdf"))
grid.arrange(pPca, pSne, ncol = 1, top = population)
dev.off()