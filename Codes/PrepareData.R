source('Codes/Functions.R')

Initialize()
set.seed(123) 

### loading data

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = F)
ColumnsToDrop <- c('id','X','X.1','X.2','ic','ev','seq','DB','annotation','rnaType')
designMat <- TotalMatrixWithStruct[,!colnames(TotalMatrixWithStruct) %in% ColumnsToDrop ]

### making all the designMatrices 

designMat_minimal <- subset(designMat, select=c(chr:strand,label))
designMat_minimal_SecondStruct <- subset(designMat, select = (label:G111))
designMat_SecondStruct <- subset(designMat, select= -(AAAA:UUUU))
designMat_NoSecStr <- subset(designMat, select=-(FreeEnergy:G111))

### spliting designMatrices to label and features

DesignMetrices <- list(designMat_minimal, designMat_minimal_SecondStruct, designMat_SecondStruct, designMat, designMat_NoSecStr)
SplitedDesignMetrices <- lapply(DesignMetrices, SplitLabelFromFeatures)
Labels <- lapply(SplitedDesignMetrices, getLabel)
Features <- lapply(SplitedDesignMetrices, getFeatures)


### Spliting designMat to train and test

SplitList <- lapply(Labels, function(aListOflables) sample.split(aListOflables , SplitRatio = 0.75) )

trainSet <- sapply(1:length(DesignMetrices),
                   function(i){ 
                     atrainSet = subset(DesignMetrices[[i]], SplitList[[i]])
                     colsToBeFactorized = c('chr', 'strand', 'label')
                     colsToBeFactorized = colsToBeFactorized[colsToBeFactorized %in% colnames(atrainSet)]
                     if(length(colsToBeFactorized) ) atrainSet[colsToBeFactorized] = lapply(atrainSet[colsToBeFactorized], factor)
                     atrainSet} , simplify = F)


testSet <- sapply(1:length(DesignMetrices),
                  function(i){ 
                    atestSet = subset(DesignMetrices[[i]], !SplitList[[i]])
                    colsToBeFactorized = c('chr', 'strand', 'label')
                    colsToBeFactorized = colsToBeFactorized[colsToBeFactorized %in% colnames(atestSet)]
                    if(length(colsToBeFactorized) ) atestSet[colsToBeFactorized] = lapply(atestSet[colsToBeFactorized], factor)
                    atestSet} , simplify = F)

