suppressPackageStartupMessages(library(randomForest))
set.seed(101)

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'trained_mod', 't', 1, 'character',
    'output', 'o', 1, 'character',
    'help', 'h', 0, 'logical'
  ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  } else if (is.null(opt$output)){
    opt <- list(ARGS=NULL)
  }
} else {
  opt <- list(ARGS=NULL)
}
# load model
load(opts$trained_mod)

# read/process feature data frame
featDf <- read.delim(opt$input)
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
rf.model <- randomForest(class ~ . - gene - strata, 
                         strata=featDf$strata, sampsize=geneCounts,
                         data=featDf)
result <- predict(rf.model, type="prob")
finalDf <- cbind(featDf, result)
write.table(finalDf, "../output/gene_features_1_26_2017/oob_predictions.txt", sep='\t')
