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

COUNTS_DATA_DIR = 'Data/smRNA_counts'
filesPath = list.files(COUNTS_DATA_DIR, full.names = T)
countsFiles = lapply(filesPath, read.delim, header=T)
names(countsFiles) <- substr(list.files(COUNTS_DATA_DIR),1,2)

countsFiles <- lapply(countsFiles, function(aFile) {
  rownames(aFile) = aFile$X
  aFile$X = NULL
  aFile
})



## SKBR==SKBR3  & CNLM1A==LM1A
### HCC1395 needs to be removed 

RefineCellLineNames <-  function(Files) {
  
  EVcellLines <- toupper(gsub('\\..*','',colnames(Files[['EV']])))
  ICcellLines <- toupper(gsub('\\..*','',gsub('\\_.*','',colnames(Files[['IC']]))))
  
  ICcellLines[ICcellLines=='SKBR'] <- 'SKBR3'
  EVcellLines[EVcellLines=='CNLM1A'] <- 'LM1A' 
  
  colnames(Files[['EV']]) <- EVcellLines
  colnames(Files[['CM']]) <- EVcellLines
  colnames(Files[['IC']]) <- ICcellLines
  
  Files <- lapply(Files, function(x) x[,!colnames(x) %in% c('HCC1395','HCC1395.1')] )
  Files
}

Counts <- RefineCellLineNames(countsFiles)



### filtering Counts based on Valid IDs defined by Ali
ValidIDs <- designMat$id
Counts <- lapply(Counts, function(x) x[rownames(x) %in% ValidIDs,])



#### Filtering noisy (false) regions (samples) based on IC expression value ####

IC_NOISE_COUNTS_THRESHOLD <- quantile(rowSums(Counts[['IC']]), 0.75)
IC_FILTER_PASSED_INDICES <- rowSums(Counts[['IC']]) > IC_NOISE_COUNTS_THRESHOLD

print(paste0('Percentage of Passed Regions: %' ,
       round(sum(IC_FILTER_PASSED_INDICES)/length(IC_FILTER_PASSED_INDICES),3)*100))

FilterCounts <- lapply(Counts, function(x) x[IC_FILTER_PASSED_INDICES,])
names(FilterCounts) <- names(Counts)

EV_EXPORTED_THRESHOLD = 20
EV_FILTER_PASSED_INDICES <- rowSums(FilterCounts[['EV']]) > EV_EXPORTED_THRESHOLD


LabelMat <- data.frame(ic=rowSums(FilterCounts[['IC']]), ev=rowSums(FilterCounts[['EV']]), 
           label=ifelse(EV_FILTER_PASSED_INDICES,'YES','NO'))

write.csv(LabelMat, file = 'Data/TertiaryLabelsMat.csv', quote = F, na = 'NA', row.names = T, col.names = T)





##### making the final design matrix based on new IC_threshold ####

designMatNgram <- read.csv('Data/NgramLen4DesignMatrix.csv',stringsAsFactors = F)
FilterDesignMat <- designMatNgram[IC_FILTER_PASSED_INDICES,]
write.csv(FilterDesignMat, 'Data/filterNgramDesignMat.csv',quote = F, na = 'NA', row.names = T, col.names = T)

## checking if filtered regions are compatible between DesignMatrix and label-Matrix
sum(rownames(FilterCounts[[1]])!=FilterDesignMat$X)
sum(LabelMat$X!=FilterDesignMat$X)





## checking length distribution 
## filtering regions with unusual lengths

bed <- read.table('Data/Annotation/combined.cluster.scoref10.bed')
colnames(bed) <- c('chr', 'start', 'end', 'id', 'score','strand')
bed$width <- bed$end - bed$start

LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)
FilteredBed <- merge(LabelMat,bed,by.x='X',by.y='id',all.x=T)


lenDisChr <- ggplot(bed, aes(x=width,y=chr))+geom_boxplot(aes(fill=chr))+theme(legend.position = "none")+ggtitle('width distribution(raw bed)')
lenDisChrFilt <- ggplot(FilteredBed, aes(x=width,y=chr))+geom_boxplot(aes(fill=chr))+theme(legend.position = "none")+ggtitle('width distribution(filtered bed)')


DrawBarPlot <- function(LabelColumn){
  ggplot(data.frame(LabelColumn), aes(x=Var1,y=Freq,color='black'))+
    geom_bar(stat = 'identity',color="dark blue", fill="cadetblue2",width=0.4)+xlab('label')+theme_bw()
}

DrawBoxPlot <- function(DataTable, Title){
  ggplot(DataTable, aes(x=label,y=width))+
    geom_boxplot(aes(fill=label))+theme_bw()+
    ggtitle(Title)
}

DrawFeatureDistribution_LabelBased(subset(TotalMatrixWithStruct, select=c(a, u, c, g, label)))



LENGTH_MIN = 18
LENGTH_MAX = 500

pdf('RegionLengthDis.pdf',width=9)
grid.arrange(lenDisChr, lenDisChrFilt,nrow=1,ncol=2)
grid.arrange(
  DrawBoxPlot(FilteredBed, 'width distribution'), 
  DrawBarPlot(table(FilteredBed$label)),
  DrawBoxPlot(subset(FilteredBed, width > LENGTH_MAX), 'width-Dis(>500)'), 
  DrawBarPlot(table(FilteredBed[FilteredBed$width>LENGTH_MAX ,'label'])),
  ncol=2,nrow=2
)
grid.arrange(
  DrawBoxPlot(subset(FilteredBed, width < LENGTH_MIN ), 'width-Dis(<15)'), 
  DrawBarPlot(table(FilteredBed[FilteredBed$width < LENGTH_MIN,'label'])),
  DrawBoxPlot(subset(FilteredBed, width < LENGTH_MAX & width > LENGTH_MIN), 'width-Dis(15:500)'), 
  DrawBarPlot(table(FilteredBed[FilteredBed$width<LENGTH_MAX &  FilteredBed$width>LENGTH_MIN,'label'])),
  ncol=2,nrow=2
)
dev.off()


### down-limit leads to loosing many labels 
print (paste0('number of labels without min-limit: ',
              nrow(subset(FilteredBed, width<LENGTH_MAX))))

print(paste0('number of labels with min-limit: ', 
             nrow(subset(FilteredBed, width<LENGTH_MAX & width>LENGTH_MIN))))

LENGTH_MIN = 18


#### making the final feature+label matrix and filtering it based on region length

designMat <- read.csv('Data/filterNgramDesignMat.csv',stringsAsFactors = F)
LabelMat <- read.csv('Data/TertiaryLabelsMat.csv', stringsAsFactors = F)
TotalMatrix <- merge(designMat, LabelMat,'X','X',all.x=T,all.y = T)
TotalMatrix <- TotalMatrix[FilteredBed$width < LENGTH_MAX & FilteredBed$width > LENGTH_MIN ,]
write.csv(TotalMatrix, 'Data/MergedDesignMatLabel_LenFilter.csv',quote = F, na = 'NA', row.names = F, col.names = T)





#TotalMatrixWithStruct = read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv',stringsAsFactors = F)


### visualize feature distribution 
pdf('plots/FeatureCor.pdf',width = 10, height = 10)
pheatmap(cor(subset(designMat,select=c(a,c,u,g))))
pheatmap(cor(Nucleotide_KmerMatrix),show_rownames = F, show_colnames = F)
pheatmap(cor(RNAsecondaryStructFeatures[,6:ncol(RNAsecondaryStructFeatures)-1]),main = 'DotBracket-Feature-correlation')
dev.off()


### label-based distributions
DrawFeatureDistribution_LabelBased(subset(TotalMatrixWithStruct, select=c(a, u, c, g, label)))
DrawFeatureDistribution_LabelBased(subset(TotalMatrixWithStruct, select=c(X0, X1, label)))
DrawFeatureDistribution_LabelBased(subset(TotalMatrixWithStruct, select=c(All_Possible_Kmers,'label')))
DrawFeatureDistribution_LabelBased(subset(TotalMatrixWithStruct, 
                                          select=c(colnames(TotalMatrixWithStruct)[278:ncol(TotalMatrixWithStruct)-1],'label')))


pdf('plots/DotBracketTriMerDistr_labelBasedPerc.pdf',height = 12)
pdf('plots/NucKmerDisLabelBasedPerc.pdf',height = 35, width = 6)
dev.off()


