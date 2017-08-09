library(randomForest)

# read/process feature data frame
#featDf <- read.delim("../output/mc3/features.txt")
featDf <- read.delim("../output/gene_features_1_26_2017/features_merged.txt")
#featDf <- read.delim("../output/gene_features_1_26_2017/snvbox/features.txt")
row.names(featDf) <- featDf$UID
featDf <- featDf[!duplicated(featDf$ID),]
colSelect <- -which(names(featDf) %in% c("UID", "ID"))
featDf <- featDf[, colSelect]

# set up gene column for use as strata in random forest
geneVec <- as.character(featDf$gene)
featDf$gene <- geneVec
featDf$strata <- geneVec
featDf[featDf["class"]=="passenger", "strata"] <- "not_driver"
featDf$strata<- factor(featDf$strata, levels=unique(featDf$strata))
geneCounts <- table(featDf$strata)

# normalize high counts
#notCt <- geneCounts["not_driver"]
medCount <- median(geneCounts)
geneCounts[geneCounts > medCount] <- medCount
totDriverCts <- sum(geneCounts[-which(names(geneCounts) %in% c("not_driver"))])
geneCounts["not_driver"] <- totDriverCts

# fill NA values
mycols <- names(featDf)
featureCols <- mycols[2:(length(mycols)-1)]
for(i in featureCols){
  featDf[is.na(featDf[,i]), i] <- mean(featDf[,i], na.rm=TRUE)
}

# perform random forest
#rf.model <- randomForest(class ~ . - gene - strata, 
                         #strata=featDf$strata, sampsize=geneCounts,
                         #sampsize=c(6000, 6000),
                         #data=featDf)
rf.model <- randomForest(class ~ . - gene, 
                         #strata=featDf$strata, sampsize=geneCounts,
                         #sampsize=c(6000, 6000),
                         data=featDf)
result <- predict(rf.model, type="prob")
finalDf <- cbind(featDf, result)
write.table(finalDf, "../output/gene_features_1_26_2017/oob_predictions.txt", sep='\t')

# get the list of genes
driverGenes <- sample(unique(featDf[featDf$class=="driver",]$gene))
driverLen <- length(driverGenes)
passengerGenes <- sample(unique(featDf[featDf$class=="passenger",]$gene))
passengerLen <- length(passengerGenes)

# add columns for driver/passenger score
featDf[,c('driver', 'passenger')] <- 0

# handle indexing
startIxD <- 1
startIxP <- 1
nfold <- 10
for (i in 1:nfold){
  # get list of genes to use for drivers
  endIxD <- i * driverLen / nfold 
  testIxD <- startIxD:endIxD
  testDGenes <- as.character(driverGenes[testIxD])
  trainDGenes <- as.character(driverGenes[-testIxD])
  startIxD <- endIxD + 1

  # get list of genes to use for passengers
  endIxP <- i * passengerLen / nfold 
  testIxP <- startIxP:endIxP
  testPGenes <- as.character(passengerGenes[testIxP])
  trainPGenes <- as.character(passengerGenes[-testIxP])
  startIxP <- endIxP + 1

  # combine list of genes
  trainGenes <- c(trainPGenes, trainDGenes)
  testGenes <- c(testPGenes, testDGenes)

  # get train and test
  trainDf <- featDf[featDf$gene %in% trainGenes,]
  testDf <- featDf[featDf$gene %in% testGenes,]
  # make sure strata has correct factor levels for training
  #tmp <- as.character(trainDf$strata)
  #trainDf$strata <- factor(tmp, levels=c(trainDGenes, 'not_driver'))

  # perform training
  rf.model <- randomForest(class ~ . - gene - driver - passenger, 
                           #strata=trainDf$strata, sampsize=geneCounts[c(trainDGenes, 'not_driver')],
                           #sampsize=c(6000, 6000),
                           data=trainDf)

  # perform test
  result <- predict(rf.model, testDf, type="prob")
  featDf[rownames(result), c('driver', 'passenger')] <- result

}
