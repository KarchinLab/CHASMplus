library(randomForest)

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'seed', 's', 1, 'integer',
    'output', 'o', 1, 'character',
    'trained_mod_file', 't', 1, 'character',
    'varimp', 'v', 1, 'character',
    'help', 'h', 0, 'logical'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  } else if (is.null(opt$input) | is.null(opt$varimp)){
    opt <- list(ARGS=NULL)
  }
} else {
  opt <- list(ARGS=NULL)
}
set.seed(opt$seed)

# read in features
featDf <- read.delim(opt$input)
row.names(featDf) <- featDf$UID
featDf <- featDf[!duplicated(featDf$ID),]
colSelect <- -which(names(featDf) %in% c("UID", "ID"))
idCol <- featDf$ID
uidCol <- featDf$UID
featDf <- featDf[, colSelect]

# setup driver pred cols
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
medCount <- median(geneCounts[-which(names(geneCounts) %in% c("not_driver"))])
geneCounts[(geneCounts > medCount) & !(names(geneCounts) %in% c("not_driver"))] <- medCount
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
passengerGenes <- sample(unique(featDf[featDf$class=="passenger",]$gene))
passengerLen <- length(passengerGenes)

# add columns for driver/passenger score
featDf[,allCols] <- 0

# train random forest
rf.model <- randomForest(class ~ . - gene - strata - driver - passenger, 
                         strata=featDf$strata, sampsize=geneCounts,
                         #sampsize=c(totDriverCts, totDriverCts),
                         data=featDf)

# save result
result <- predict(rf.model, type="prob")
#result <- rf.model$predicted
featDf[rownames(result), c('driver', 'passenger')] <- result
featDf["ID"] <- idCol
featDf["UID"] <- uidCol
write.table(featDf, opt$output, sep='\t')
write.table(rf.model$importance, opt$varimp, sep='\t')
save(rf.model, file=opt$trained_mod_file)
