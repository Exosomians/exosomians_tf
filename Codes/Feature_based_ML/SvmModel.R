source('Codes/Functions.R')

Initialize()
set.seed(123) 

trainSet <- readRDS('Data/trainSet.rds')
testSet <- readRDS('Data/testSet.rds')

BalWeight = c("NO" = 1, "YES" = 10)
unBalWeight = c("NO" = 1, "YES" = 1)

### Training Sets
trainSet_min <- trainSet[[1]]
trainSet_minSecStr <- trainSet[[2]]
trainSet_SecStr <- trainSet[[3]]
trainSet_total <- trainSet[[4]]
trainSet_noSecStr <- trainSet[[5]] 
trainSet_minKmer <- subset(trainSet_total, select= AAAA:UUUU)
trainSet_length <- subset(trainSet_total, select= c(length,label))


### Train SVM models
SVMlin_Bal <- svmModel(trainSet_total, 'linear', BalWeight)
SVMrad_Bal <- svmModel(trainSet_total, 'radial', BalWeight)
SVMlin_NoBal <- svmModel(trainSet_total, 'linear', unBalWeight)
SVMlin_NoStr <- svmModel(trainSet_noSecStr, 'linear', BalWeight)
SVMrad_NoStr <- svmModel(trainSet_noSecStr, 'radial', BalWeight)
SVMlin_SecStr_Bal <- svmModel(trainSet_SecStr, 'linear', BalWeight)
SVMlin_min_Bal <- svmModel(trainSet_min,'linear',BalWeight)
SVMlin_minSecStruct_Bal <- svmModel(trainSet_minSecStr,'linear',BalWeight)
SVMlin_minKmer_Bal <- svmModel(trainSet_minKmer, 'linear', BalWeight)
SVMlin_len_Bal <- svmModel(trainSet_length, 'linear', BalWeight)


listOfsvmModels <- list(SVMlin_Bal, SVMrad_Bal, SVMlin_NoBal, SVMlin_NoStr , SVMrad_NoStr,
                  SVMlin_SecStr_Bal, SVMlin_min_Bal, SVMlin_minSecStruct_Bal)

#saveRDS(listOfsvmModels, 'models/DNA2RNAdesignMat/listOfModels_SVM.rds')
listOfsvmModels <- readRDS('models/DNA2RNAdesignMat/listOfModels_SVM.rds')
SVMlin_minKmer_Bal <- readRDS('models/SVMlin_minKmer_Bal.rds')


### Test Sets
testSet_min <- testSet[[1]]
testSet_minSecStr <- testSet[[2]]
testSet_SecStr <- testSet[[3]]
testSet_total <- testSet[[4]]
testSet_noSecStr <- testSet[[5]]
testSet_minKmer <- subset(testSet_total, select= AAAA:UUUU)
testSet_length <- subset(testSet_total, select= length)


listOfTestSets <- list(testSet_total, testSet_total, testSet_total, testSet_noSecStr, 
                       testSet_noSecStr, testSet_SecStr, testSet_min, testSet_minSecStr )

listOfSvmPreds <- sapply(1:length(listOfTestSets), 
                      function(i){
                        predict( listOfsvmModels[[i]],  listOfTestSets[[i]])
                      }, simplify = F)


SVMlin_minKmerBalPred <- predict( SVMlin_minKmer_Bal, testSet_minKmer)
SVMlin_SVMlinlenBalPred <- predict( SVMlin_len_Bal, testSet_length)


draw_confusion_matrix(confusionMatrix(SVMlin_minKmerBalPred , listOfTestSets[[1]]$label), 
                      'SVM-weighted(1:10)-Min-kmer-linearKernel', 'No','Yes')
draw_confusion_matrix(confusionMatrix(SVMlin_SVMlinlenBalPred , listOfTestSets[[1]]$label), 
                      'SVM-weighted(1:10)-length-linearKernel', 'No','Yes')



saveRDS(listOfSvmPreds, 'models/listOfSvmPreds.rds')
listOfTitles <- c('SVM-weighted(1:10)-linearKernel', 
                  'SVM-weighted(1:10)-RadialKernel', 
                  'SVM-unWeighted-linearKernel', 
                  'SVM-weighted(1:10)-noStructFeature-linearKernel', 
                  'SVM-weighted(1:10)-noStructFeature-radialKernel',
                  'SVM-weighted(1:10)-SecondStruct(NO-Kmer)-linearKernel',
                  'SVM-weighted(1:10)-Minimal-linearKernel',
                  'SVM-weighted(1:10)-Min-SecStruct-linearKernel')

#### SVM Model Evaluation

pdf('plots/SVMresults_DNA2RNA.pdf')
sapply(1:length(listOfSvmPreds), 
       function(i){
         draw_confusion_matrix(
           confusionMatrix(
             data = listOfSvmPreds[[i]] , 
             reference = listOfTestSets[[i]]$label), 
           
           listOfTitles[i], 'No','Yes')} )

dev.off()


trainSet_minKmer <- subset(trainSet_total, select= c(AAAA:UUUU, label))
SVMlin_minKmer_Bal <- svmModel(trainSet_minKmer, 'linear', BalWeight)
saveRDS(SVMlin_minKmer_Bal, 'models/SVMlin_minKmer_Bal.rds')

predict( listOfsvmModels[[i]],  listOfTestSets[[i]])
testSet_minKmer <- subset(testSet_total, select= AAAA:UUUU)



