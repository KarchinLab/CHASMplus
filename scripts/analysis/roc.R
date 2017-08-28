suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(dplyr))

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'chasm2', 'c', 1, 'character',
    'poslabel', 'p', 1, 'character',
    'output', 'o', 1, 'character',
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


benchmark_df <- read.delim(opt$input, stringsAsFactors=F)
# Flip scores which are more driver-like with smaller scores
benchmark_df['1-CHASM'] <- 1-benchmark_df['CHASM']
benchmark_df['1-SIFT'] <- 1-benchmark_df['SIFT']
benchmark_df['FATHMM'] <- -benchmark_df['FATHMM']
# merge in chasm2 genome scores
chasm2_df <- read.delim(opt$chasm2, stringsAsFactors=F)
chasm2_df <- chasm2_df %>% 
              select(gene, driver.score) %>% 
              distinct()
benchmark_df <- left_join(benchmark_df, chasm2_df, by=c('Hugo_Symbol'='gene'))
# add CHASM2 genome score
benchmark_df['CHASM2_genome'] <- benchmark_df['CHASM2'] * benchmark_df['driver.score']
# establish label for positive class
benchmark_df['y'] <- as.integer(benchmark_df['class']==opt$poslabel)

# prediction methods
pred_methods <- c('VEST', 'CADD',  'FATHMM', 
                  '1-SIFT', 'MutationAssessor', 'REVEL', 'MCAP',
                  'ParsSNP', '1-CHASM',  
                  'Polyphen2_hdiv', 'Polyphen2_hvar', 
                  'CanDrA', 'CanDrA.plus',
                  'CHASM2', 'CHASM2_genome')
# calculate the AUC's
aucs <- c()
for (m in pred_methods) {
  tmp_data <- benchmark_df[!is.na(benchmark_df[,m]),]
  tmp_auc <- as.numeric(auc(tmp_data$y, tmp_data[,m]))
  aucs <- c(aucs, tmp_auc)
}
result_df <- data.frame(method=pred_methods, auc=aucs)

# figure out which is the best auc within different variants of the same method
# polyphen
polyphen <- c('Polyphen2_hvar', 'Polyphen2_hdiv')
idx_max <- which.max(result_df[result_df$method %in% polyphen, 'auc'])
best_polyphen <- polyphen[idx_max]
# candra
candra <- c('CanDrA', 'CanDrA.plus')
idx_max <- which.max(result_df[result_df$method %in% candra, 'auc'])
best_candra <- candra[idx_max]
# CHASM2
chasm2 <- c('CHASM2', 'CHASM2_genome')
idx_max <- which.max(result_df[result_df$method %in% chasm2, 'auc'])
best_chasm2 <- chasm2[idx_max]

# test for statistical significance
method_compare_vec <- c(pred_methods[1:(length(pred_methods)-6)], c(best_polyphen, best_candra, best_chasm2))
pvals <- c()
for (m in method_compare_vec) {
  tmp_data <- benchmark_df[!is.na(benchmark_df[,m]) & !is.na(benchmark_df[,best_chasm2]),]
  tmp_auc_other <- auc(tmp_data$y, tmp_data[,m])
  tmp_auc_chasm2 <- auc(tmp_data$y, tmp_data[,best_chasm2])
  roc_test <- roc.test(tmp_auc_chasm2, tmp_auc_other, method='delong')
  pvals <- c(pvals, roc_test$p.value)
}
result_df[,'pvalue'] <- NA
row.names(result_df) <- result_df$method
result_df[method_compare_vec, 'pvalue'] <- pvals

# save results
write.table(result_df, opt$output, sep='\t')
