# This script will do the hierarchical clustering and using dyamic treecut to estimate the optimal number of clusters
#clear workspace
rm(list = ls(all = TRUE))

# user define the following parameters
session <- c('rfMRI_REST2_LR')

baseResults <- sprintf("/home/data/Projects/HCP/results/%s", session)
baseFig <- sprintf("/home/data/Projects/HCP/figs/%s", session)
#winFile <- sprintf("zfeatureFC_winFullCor_%s.mat", session)
#varName <- 'zfeatureWin'
# normFeature1_winFullCor_session1.mat normFeature1
winFile <- sprintf("norm2FeatureWin_%s.mat", session)
varName <- 'norm2FeatureWin'


## The followin section dosen't need to change ##


# Load matlab library
library(R.matlab)

# Read in matlab data
fileName <- file.path(baseResults, winFile)
matData <- readMat(fileName)

# Copy over matrix
win <- matData[[varName]]

# Check dimensions
dim(win)

# import daynamic tree cut library
library(dynamicTreeCut)

# clustering and cut the tree for windows of all seeds 
dist <- dist(win, method='euclidean')
dendro <- hclust(dist,method='ward')
clustIndx <- cutreeDynamic(dendro, minClusterSize=2, 	distM=as.matrix(dist))
		
# save the index file
indxFileName <- file.path(baseResults,sprintf("clustIndx_%s.txt", winFile[1:end-4]))
write.table(clustIndx,file=indxFileName,row.names=FALSE,col.names=FALSE,qmethod="double")
	
# plot and save the dendrogram plot
figName <- file.path(baseFig,sprintf("dendro_%s.png", winFile[1:end-4]))
png(filename = figName , width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
figTitle <- sprintf("denfro_%s", winFile)
plot(dendro,main=figTitle,ylab="Heights",xlab='feature windows', labels=F)
dev.off()
		
		

		
		
							
		

