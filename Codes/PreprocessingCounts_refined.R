# In this file, we will read the data,
#   Try to filter the noisy (false) regions (samples/sequences) based on their IC expression value,
#   Then, label the data, i.e. which sample/sequence is exported via exosome
#   Then, check the width distribution of the RNA regions
#   Then, filter the label and design matrix based on width distributions

#   Input: Raw data files (smRNA_counts) 
#   Output: Labels Matrix


source('Codes/Functions.R')
Initialize()


#### Reading raw data from file ####

COUNTS_DATA_DIR = 'Data/oldDataRefined/'
filesPath = list.files(COUNTS_DATA_DIR,pattern = '*.txt' ,full.names = T)
countsFiles = lapply(filesPath, read.delim, header=F)
names(countsFiles) <- substr(list.files(COUNTS_DATA_DIR,pattern = '*.txt' ),1,2)

designMat = read.csv('Data/oldDataRefined/1_PrimaryDesignMat.csv')

countsFiles <- lapply(countsFiles, function(aFile) {
  aFile <- subset(aFile, select=c(V1, V2, V3, V4, V5))
  colnames(aFile) = c('chr', 'start', 'end', 'coverage', 'length')
  aFile$id <- paste0(aFile$chr, '_', aFile$start, '_', aFile$end)
  aFile <- subset(aFile, id %in% designMat$id)
})

## checjing if the regions are in order
lapply(countsFiles, function(x) sum(x$id != designMat$id))
designMat$ic = countsFiles[['IC']]$coverage
designMat$ev = countsFiles[['EV']]$coverage


## defining labels
designMat$label <- ifelse( designMat$ev>quantile(designMat$ev, 0.75), 'YES', 'NO')
View(head(designMat,10))
write.csv(designMat, 'Data/oldDataRefined/2_PrimaryDesignMat_label.csv')

check = read.csv('Data/oldDataRefined/2_PrimaryDesignMat_label.csv')
head(check)
