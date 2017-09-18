suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(dplyr))

# get command line args
if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'input', 'i', 1, 'character',
    'chasm2', 'c', 1, 'character',
    'list_cgc', 'l', 2, 'character',
    'mutations', 'm', 2, 'character',
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

# function to read CGC genes
read_cgc <- function(path){
  # read data
  cgc_df <- read.delim(path, stringsAsFactors=F) 

  # keep only somatic missense genes
  cgc_genes <- cgc_df  %>%
    filter(!grepl('Mis', Mutation.Types) & Tumour.Types.Somatic.!='') %>%
    select(Gene.Symbol) %>%
    distinct %>%
    as.vector

  return(cgc_genes[,'Gene.Symbol'])
}

# merge the recurrence of a mutation into the 
merge_recurrence <- function(bench_df, mut_df) {
  # count the recurence of the mutation
  recur <- mut_df %>%
              filter(Variant_Classification=='Missense_Mutation') %>%
              group_by(Hugo_Symbol, HGVSp_Short) %>%
              summarise(recurrence=n())

  # merge info back into benchmark
  bench_df <- bench_df %>%
                left_join(recur, by=c('Hugo_Symbol', 'HGVSp_Short')) %>%
                as.data.frame

  return(bench_df)
}



benchmark_df <- read.delim(opt$input, stringsAsFactors=F)
# Flip scores which are more driver-like with smaller scores
benchmark_df['1-CHASM'] <- 1-benchmark_df['CHASM']
benchmark_df['1-SIFT'] <- 1-benchmark_df['SIFT']
benchmark_df['FATHMM'] <- -benchmark_df['FATHMM']

poslabel <- opt$poslabel
if (poslabel!='mc3') {
  # merge in chasm2 genome scores
  chasm2_df <- read.delim(opt$chasm2, stringsAsFactors=F)
  chasm2_df <- chasm2_df %>% 
                select(gene, driver.score) %>% 
                distinct()
  benchmark_df <- left_join(benchmark_df, chasm2_df, by=c('Hugo_Symbol'='gene'))
  # add CHASM2 genome score
  benchmark_df['CHASM2_genome'] <- benchmark_df['CHASM2'] * benchmark_df['driver.score']
}


# establish label for positive class
if (poslabel=='iarc_tp53') {
  benchmark_df['y'] <- as.integer(benchmark_df['class']<50)
} else if (poslabel=='activating') {
  benchmark_df <- benchmark_df[benchmark_df$class %in% c('activating', 'neutral'),]
  benchmark_df['y'] <- as.integer(benchmark_df['class']==opt$poslabel)
} else if (poslabel=='mc3') {
  # merge in chasm2 scores
  chasm2_df <- read.delim(opt$chasm2, stringsAsFactors=F)
  chasm2_df <- chasm2_df %>%
                select(UID, CHASM2, CHASM2_genome, 
                       CHASM2_pval, CHASM2_genome_pval,
                       CHASM2_qval, CHASM2_genome_qval)
  benchmark_df <- benchmark_df %>% 
                    left_join(chasm2_df, by=c('UID'))

  # define the class labels
  cgc_genes <- read_cgc(opt$list_cgc)
  mut_df <- read.delim(opt$mutations, stringsAsFactors=F)
  benchmark_df <- merge_recurrence(benchmark_df, mut_df)
  is_cgc <- benchmark_df[,'Hugo_Symbol'] %in% cgc_genes
  is_recur <- benchmark_df['recurrence'] > 1
  benchmark_df['y'] <- as.integer(is_cgc & is_recur)
} else if (poslabel=='msk_impact') {
  benchmark_df['y'] <- as.integer(benchmark_df$class %in% c('Oncogenic', 'Likely Oncogenic'))
} else {
  benchmark_df['y'] <- as.integer(benchmark_df['class']==opt$poslabel)
}

# prediction methods
print('Calculating AUC . . .')
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
print('Finished calculating AUC.')

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
  print(paste('Comparing', m))
  tmp_data <- benchmark_df[!is.na(benchmark_df[,m]) & !is.na(benchmark_df[,best_chasm2]),]
  tmp_auc_other <- auc(tmp_data$y, tmp_data[,m])
  tmp_auc_chasm2 <- auc(tmp_data$y, tmp_data[,best_chasm2])
  roc_test <- roc.test(tmp_auc_chasm2, tmp_auc_other, method='delong')
  pvals <- c(pvals, roc_test$p.value)
  print('finished')
}
result_df[,'pvalue'] <- NA
row.names(result_df) <- result_df$method
result_df[method_compare_vec, 'pvalue'] <- pvals

# save results
write.table(result_df, opt$output, sep='\t')
