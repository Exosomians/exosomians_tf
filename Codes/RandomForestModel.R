source('Codes/Functions.R')
source('Codes/PrepareData.R')

library("pROC")
library('ROCR')

Initialize()
set.seed(123) 
h2o.init(ip = 'localhost', port = 54321, nthreads= detectCores()-4, max_mem_size = '48g')

designMatNames <- c('designMat_minimal', 'designMat_minimal_SecondStruct', 
                    'designMat_SecondStruct', 'designMat', 'designMat_NoSecStr')

modelNames <- gsub('designMat','RF',designMatNames)
featureNames <- lapply(Features, colnames)
label = 'label'

trainSet.h2o <- lapply(trainSet, as.h2o )
testSet.h2o <- lapply(testSet, function(x) as.h2o(subset(x, select=-label)))

#### train Random Forest

listOfRFModels <- sapply(1:length(trainSet), 
       function(i){ 
         h2o.randomForest(
           x = featureNames[[i]] , 
           y = label, 
           training_frame = trainSet.h2o[[i]], 
           balance_classes = T, 
           ntrees = 100, 
           max_depth = 20, 
           seed = 123)}, 
       simplify = F)

names(listOfRFModels) <- modelNames

#### RF performance

listOfFeatureImportance <- lapply(listOfRFModels, h2o.varimp )
names(listOfFeatureImportance) <- modelNames
lapply(listOfFeatureImportance, function(x) x[x$scaled_importance>0.2,])

listOfRFPreds <- sapply(1:length(listOfRFModels), 
                        function(i) as.data.frame( h2o.predict(listOfRFModels[[i]], testSet.h2o[[i]]) ), 
                        simplify = F)


pdf('plots/RFresults_DNA2RNA.pdf')
sapply(1:length(listOfRFModels), 
       function(i){
         draw_confusion_matrix(
           confusionMatrix(
             data=listOfRFPreds[[i]]$predict, 
             reference = testSet[[i]]$label),
           modelNames[i], 'No','Yes', NULL,'rosybrown1')
         
       } )
dev.off()


ListOfRF_ROCs <- sapply(1:length(listOfRFModels), 
                        function(i){
                          ROC = roc(
                            predictor = MakeNumeric(listOfRFPreds[[i]]$predict),
                            response = MakeNumeric(testSet[[i]]$label))
                          return(ROC$auc) }
                        , simplify = F)


## what?? 
ListOfRF_AUCs <- sapply(1:length(listOfRFModels), 
                        function(i){
                          AUC = ROCR::prediction(
                            predictions = MakeNumeric(listOfRFPreds[[i]]$predict),
                            labels = MakeNumeric(testSet[[i]]$label))
                          return(ROCR::performance(AUC,'auc')) }
                        , simplify = F)


