suppressPackageStartupMessages(library(randomForest))

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'seed', 's', 1, 'integer',
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

# permute class labels
allgenes <- sample(unique(as.character(featDf$gene)))
gene_cts <- table(as.character(featDf$gene))
dgenes <- unique(as.character(featDf[featDf["class"]=="driver",]$gene))
featDf_subset <- featDf[featDf$gene %in% dgenes,]
mytable <- table(as.character(featDf_subset$gene), as.character(featDf_subset$class))
driver_frac <- mytable[,"driver"] / (mytable[,"driver"] + mytable[,"passenger"])
tot_num_driver <- sum(mytable[,"driver"])
featDf["class"] <- "passenger"  # reset everything to passengers
num_driver_sum <- 0

i <- 1
j <- 1
while (num_driver_sum < tot_num_driver) {
  # start over on driver fractions once
  # passed over all that were observed
  if (i>length(driver_frac)){
    i <- 1 
  }
  # randomly assign the same fraction of drivers
  # as actually used within the observed data
  mygene <- allgenes[j]
  tot_mut <- gene_cts[mygene]
  num_driver <- ceiling(tot_mut * driver_frac[i])
  driver_vec <- sample(c(rep("driver", num_driver), rep("passenger", tot_mut-num_driver)))
  featDf[featDf["gene"]==mygene,"class"] <- driver_vec

  # update counters
  i <- i+1
  j <- j+1
  # update number of assigned drivers
  num_driver_sum <- num_driver_sum + num_driver
}

#featDf[,"class"] <- sample(featDf$class)

# fix factor
featDf[,"class"] <- as.factor(featDf$class)

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
                         data=featDf)

# save result
#result <- predict(rf.model, featDf, type="prob")
#featDf[rownames(result), c('driver', 'passenger')] <- result
#featDf["ID"] <- idCol
#featDf["UID"] <- uidCol
write.table(rf.model$importance, opt$varimp, sep='\t')
