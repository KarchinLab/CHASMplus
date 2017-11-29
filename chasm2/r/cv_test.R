suppressPackageStartupMessages(library(randomForest))
set.seed(101)

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'trained_mod_dir', 't', 1, 'character',
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
  print("Getopt library is not installed")
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
driverCols <- c('driver')
allCols <- c(driverCols, 'passenger')

# fill na
featureCols <- which(sapply(featDf, is.numeric))
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

# add random other columns used in training
featDf["strata"] <- "passenger"

# handle indexing
startIxD <- 1
startIxP <- 1
nfold <- 10
for (i in 1:nfold){
  # read trained model
  #print('Reading model ...')
  myfile <- paste(paste(opt$trained_mod_dir, "train", sep="/"), i, ".Rdata", sep="")
  load(myfile)
  #print('Done reading model ...')

  # establish test data frame
  testDf <- featDf[!featDf["gene"] %in% train_genes,]

  #print("Making predictions")
  # update predictions
  result <- predict(rf.model, testDf, type="prob")
  featDf[rownames(result), c('driver', 'passenger')] <- result
  #print("Finished predictions")
}
# add id information back
featDf["ID"] <- idCol
featDf["UID"] <- uidCol

# save output
write.table(featDf, opt$output, sep='\t')
