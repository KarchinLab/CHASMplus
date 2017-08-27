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
load(opt$trained_mod)

# read/process feature data frame
featDf <- read.delim(opt$input)
row.names(featDf) <- featDf$UID
featDf <- featDf[!duplicated(featDf$ID),]
colSelect <- -which(names(featDf) %in% c("UID", "ID"))
idCol <- featDf$ID
uidCol <- featDf$UID
featDf <- featDf[, colSelect]

# set up gene column for use as strata in random forest
geneVec <- as.character(featDf$gene)
featDf$gene <- geneVec
featDf$strata <- geneVec
featDf[featDf["class"]=="passenger", "strata"] <- "not_driver"
featDf$strata<- factor(featDf$strata, levels=unique(featDf$strata))

# fill NA values
mycols <- names(featDf)
featureCols <- mycols[2:(length(mycols)-2)]
for(i in featureCols){
  featDf[is.na(featDf[,i]), i] <- mean(featDf[,i], na.rm=TRUE)
}

# perform random forest
featDf['driver'] <- 0
featDf['passenger'] <- 0
result <- predict(rf.model, featDf, type="prob")
featDf[, c('driver', 'passenger')] <- result
featDf["ID"] <- idCol
featDf["UID"] <- uidCol
write.table(featDf, opt$output, sep='\t')
