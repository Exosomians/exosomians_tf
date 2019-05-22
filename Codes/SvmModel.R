source('Codes/Functions.R')
source('Codes/PrepareData.R')

Initialize()
set.seed(123) 


BalWeight = c("NO" = 1, "YES" = 10)
unBalWeight = c("NO" = 1, "YES" = 1)

### Training Sets
trainSet_min <- trainSet[[1]]
trainSet_minSecStr <- trainSet[[2]]
trainSet_SecStr <- trainSet[[3]]
trainSet_total <- trainSet[[4]]
trainSet_noSecStr <- trainSet[[5]] 


### Train SVM models
SVMlin_Bal <- svmModel(trainSet_total, 'linear', BalWeight)
SVMrad_Bal <- svmModel(trainSet_total, 'radial', BalWeight)
SVMlin_NoBal <- svmModel(trainSet_total, 'linear', unBalWeight)
SVMlin_NoStr <- svmModel(trainSet_noSecStr, 'linear', BalWeight)
SVMrad_NoStr <- svmModel(trainSet_noSecStr, 'radial', BalWeight)
SVMlin_SecStr_Bal <- svmModel(trainSet_SecStr, 'linear', BalWeight)
SVMlin_min_Bal <- svmModel(trainSet_min,'linear',BalWeight)
SVMlin_minSecStruct_Bal <- svmModel(trainSet_minSecStr,'linear',BalWeight)


listOfModels <- list(SVMlin_Bal, SVMrad_Bal, SVMlin_NoBal, SVMlin_NoStr , SVMrad_NoStr,
                  SVMlin_SecStr_Bal, SVMlin_min_Bal, SVMlin_minSecStruct_Bal)

saveRDS(listOfModels, 'models/DNA2RNAdesignMat/listOfModels_SVM.rds')


### Test Sets
testSet_min <- testSet[[1]]
testSet_minSecStr <- testSet[[2]]
testSet_SecStr <- testSet[[3]]
testSet_total <- testSet[[4]]
testSet_noSecStr <- testSet[[5]]


listOfTestSets <- list(testSet_total, testSet_total, testSet_total, testSet_noSecStr, 
                       testSet_noSecStr, testSet_SecStr, testSet_min, testSet_minSecStr )

listOfPreds <- sapply(1:length(listOfTestSets), 
                      function(i){
                        predict( listOfModels[[i]],  listOfTestSets[[i]])
                      }, simplify = F)


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
sapply(1:length(listOfPreds), 
       function(i){
         draw_confusion_matrix(
           confusionMatrix(
             data = listOfPreds[[i]] , 
             reference = listOfTestSets[[i]]$label), 
           
           listOfTitles[i], 'No','Yes')} )

dev.off()


