# In this file, we will read the data,
#   Try to do some explorary data analysis on them,
#   Then filter the noisy (false) regions (samples/sequences) based on their IC expression value,
#   Then label the data, i.e. which sample/sequence is exported via exosome

#   Input: Raw data files (smRNA_counts) 
#   Output: Labels Matrix

# Make sure that 'libssl-dev', 'libgit2-dev', and 'libxml2-dev' libraries are installed on your OS.
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
# install.packages("tidyverse")
source('Codes/Functions.R')
Initialize()


#### Reading raw data from file ####
COUNTS_DATA_DIR = 'Data/smRNA_counts/'
filesPath = list.files(COUNTS_DATA_DIR, full.names = T)
countsFiles = lapply(filesPath, read.delim, header=T)
countsFiles = lapply(countsFiles, function(aFile) {
  rownames(aFile) = aFile$X
  aFile$X = NULL
  tempColnames = gsub(colnames(aFile), pattern = '[.].*|_.*', replacement = '')
  tempColnames = paste(tempColnames, FindRepetitionIndex(tempColnames), sep = '_')
  colnames(aFile) = tempColnames
  aFile
})

filesName = list.files(COUNTS_DATA_DIR, full.names = F)
names(countsFiles) = gsub(filesName, pattern = '_.*', replacement = '')



icCounts = countsFiles$IC

#### Getting some summaries from IC counts ####
icSummaries = CalcSummaryForCellLines(icCounts)

#### Focus on normal cell line (HUMEC) ####
icNormal = icCounts[, grepl(colnames(icCounts), pattern = 'HUMEC', ignore.case = T)]
# icNormal$HUMEC_Mean = floor(rowMeans(icNormal))
icNormal$HUMEC_Mean = round(rowMeans(icNormal), 1)


QUANTILE_LIMIT_ON_PLOT = 0.90
icNormalPlot = ggplot(icNormal) +
  geom_density(aes(x = HUMEC_Mean)) +
  xlim(range(1, quantile(icNormal$HUMEC_Mean, QUANTILE_LIMIT_ON_PLOT)))

# icNormalPlot


#### Filtering noisy (false) regions (samples) based on IC expression value ####
IC_NOISE_COUNTS_THRESHOLD = 1
icNormalFiltered = subset(icNormal, HUMEC_Mean>= IC_NOISE_COUNTS_THRESHOLD)
filterPercentage = round((1- nrow(icNormalFiltered)/nrow(icNormal)) * 100, 2)
sprintf('%0.2f%% of sequences is filtered', filterPercentage)

evCounts = countsFiles$EV
evNormal = evCounts[, grepl(colnames(evCounts), pattern = 'HUMEC', ignore.case = T)]
evNormal$HUMEC_Mean = round(rowMeans(evNormal), 1)


icPreLabeled = data.frame(id = rownames(icNormalFiltered),
                          ic = icNormalFiltered$HUMEC_Mean)
evPreLabeled = data.frame(id = rownames(evNormal),
                          ev = evNormal$HUMEC_Mean)


labelsMat = merge(icPreLabeled, evPreLabeled, all.x = T)

write.csv(labelsMat, file = 'Data/PrimaryLabelsMat.csv', quote = F, na = 'NA', row.names = F, col.names = T)

#### Labeling data (exported or not) based on EV expression value ####
EV_EXPORTED_THRESHOLD = 1
labelsMat$label = labelsMat$ev >= EV_EXPORTED_THRESHOLD
labelsMat$label[labelsMat$label==T] = 'Yes'
labelsMat$label[labelsMat$label==F] = 'No'

write.csv(labelsMat, file = 'Data/SecondaryLabelsMat.csv', quote = F, na = 'NA', row.names = F, col.names = T)
