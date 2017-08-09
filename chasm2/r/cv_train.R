suppressPackageStartupMessages(library(randomForest))
set.seed(101)

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'trained_mod_dir', 't', 1, 'character',
    'help', 'h', 0, 'logical'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  } 
} else {
  opt <- list(ARGS=NULL)
}

# read in features
featDf <- read.delim(opt$input)
row.names(featDf) <- featDf$UID
featDf <- featDf[!duplicated(featDf$ID),]
colSelect <- -which(names(featDf) %in% c("UID", "ID"))
idCol <- featDf$ID
uidCol <- featDf$UID
featDf <- featDf[, colSelect]

# setup driver pred cols
#driverCols <- c('oncogene', 'tsg')
driverCols <- c('driver')
allCols <- c(driverCols, 'passenger')

# set up gene column for use as strata in random forest
geneVec <- as.character(featDf$gene)
featDf$gene <- geneVec
featDf$strata <- geneVec
featDf[featDf["class"]=="passenger", "strata"] <- "not_driver"
featDf$strata <- factor(featDf$strata, levels=unique(featDf$strata))
geneCounts <- table(featDf$strata)

# normalize high counts
medCount <- median(geneCounts)
geneCounts[geneCounts > medCount] <- medCount
#geneCounts[] <- 1
totDriverCts <- sum(geneCounts[-which(names(geneCounts) %in% c("not_driver"))])
geneCounts["not_driver"] <- totDriverCts

# fill na
mycols <- names(featDf)
featureCols <- mycols[2:(length(mycols)-2)]
for(i in featureCols){
  featDf[is.na(featDf[,i]), i] <- mean(featDf[,i], na.rm=TRUE)
}

# get the list of genes
driverGenes <- sample(unique(featDf[featDf$class %in% driverCols,]$gene))
driverLen <- length(driverGenes)
tmp1 <- unique(featDf[featDf$class=="passenger",]$gene)
tmp2 <- setdiff(tmp1, driverGenes)
passengerGenes <- sample(tmp2)
passengerLen <- length(passengerGenes)

# add columns for driver/passenger score
featDf[,allCols] <- 0

# handle indexing
startIxD <- 1
startIxP <- 1
nfold <- 10
for (i in 1:nfold){
  # get list of genes to use for drivers
  endIxD <- ceiling(i * driverLen / nfold)
  testIxD <- startIxD:endIxD
  testDGenes <- as.character(driverGenes[testIxD])
  trainDGenes <- as.character(driverGenes[-testIxD])
  startIxD <- endIxD + 1

  # get list of genes to use for passengers
  endIxP <- ceiling(i * passengerLen / nfold)
  testIxP <- startIxP:endIxP
  testPGenes <- as.character(passengerGenes[testIxP])
  trainPGenes <- as.character(passengerGenes[-testIxP])
  # fix to remove all driver train genes
  #testPGenes <- setdiff(testPGenes, trainDGenes)
  startIxP <- endIxP + 1

  # combine list of genes
  trainGenes <- c(trainPGenes, trainDGenes)
  testGenes <- c(testPGenes, testDGenes)

  # get train and test
  trainDf <- featDf[featDf$gene %in% trainGenes,]
  testDf <- featDf[featDf$gene %in% testGenes,]
  # make sure strata has correct factor levels for training
  tmp <- as.character(trainDf$strata)
  trainDf$strata <- factor(tmp, levels=c(trainDGenes, 'not_driver'))

  # set up gene column for use as strata in random forest
  geneVec <- as.character(trainDf$gene)
  trainDf$gene <- geneVec
  trainDf$strata <- geneVec
  trainDf[trainDf["class"]=="passenger", "strata"] <- "not_driver"
  trainDf$strata <- factor(trainDf$strata, levels=unique(trainDf$strata))
  geneCounts <- table(trainDf$strata)

  # normalize high counts
  medCount <- median(geneCounts[-which(names(geneCounts) %in% c("not_driver"))])
  geneCounts[(geneCounts > medCount) & !(names(geneCounts) %in% c("not_driver"))] <- medCount
  totDriverCts <- sum(geneCounts[-which(names(geneCounts) %in% c("not_driver"))])
  geneCounts["not_driver"] <- totDriverCts

  # train random forest
  #print('training model . . .')
  rf.model <- randomForest(class ~ . - gene - strata - driver - passenger, 
                           strata=trainDf$strata, sampsize=geneCounts,
                           data=trainDf)
  #print('Done training model.')

  # save trained model
  #print('Saving model ...')
  train_genes <- unique(trainDf$gene)
  myfile <- paste(paste(opt$trained_mod_dir, "train", sep="/"), i, ".Rdata", sep="")
  save(rf.model, train_genes, file=myfile)
  #print('Done saving model ...')

  # update predictions
  result <- predict(rf.model, testDf, type="prob")
  featDf[rownames(result), c('driver', 'passenger')] <- result
}

# don't write output
#featDf["ID"] <- idCol
#featDf["UID"] <- uidCol
#write.table(featDf, opt$output, sep='\t')
#write.table(rf.model$importance, opt$varimp, sep='\t')
