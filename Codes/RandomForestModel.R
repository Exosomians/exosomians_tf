source('Codes/Functions.R')
Initialize()
h2o.init(ip = 'localhost', port = 54321, nthreads= detectCores()-4, max_mem_size = '48g')



TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = F)
LabelMat <- subset(TotalMatrixWithStruct, select=c('id','ic','ev','label'))
ColumnsToDrop <- c('id','X','X.1','ic','ev','seq','DB','annotation','rnaType')
designMat <- TotalMatrixWithStruct[,!colnames(TotalMatrixWithStruct) %in% ColumnsToDrop ]
colnames(designMat)


##### Data Preperation 
set.seed(123) 
split = sample.split(designMat$label , SplitRatio = 0.9) 
training_set = subset(designMat, split == TRUE ) 
test_set = subset(designMat, split == FALSE)


colsToBeFactorized = c('chr', 'label', 'strand')
training_set[colsToBeFactorized] = lapply(training_set[colsToBeFactorized], factor)
test_set[colsToBeFactorized] = lapply(test_set[colsToBeFactorized], factor)


features = colnames(subset(training_set,select= -label))
label = 'label'
training_set.h2o = as.h2o(training_set)
test_set.h2o = as.h2o(test_set)

####### Random Forest

RF_NoBal = h2o.randomForest(x=features, y=label, training_frame=training_set.h2o, 
                            balance_classes=F, ntrees=50, max_depth = 20, seed = 1398)

RF_Bal = h2o.randomForest(x =features, y = label, training_frame = training_set.h2o,
                       balance_classes = T, ntrees = 50, max_depth = 20, seed = 1398)

RF_Bal_200t = h2o.randomForest(x =features, y = label, training_frame = training_set.h2o,
                       balance_classes = T, ntrees = 200, max_depth = 20, seed = 1398)



library("pROC")
##############################


importance.h2o<-h2o.varimp(RF_Bal)
importance.h2o[importance.h2o$scaled_importance>0.4,]


###  RF_Bal performance
RF_NoBal_Pred <- as.data.frame(h2o.predict(RF_NoBal,test_set.h2o))
RF_Bal_Pred <- as.data.frame(h2o.predict(RF_Bal,test_set.h2o))
RF_Bal_200t_Pred <- as.data.frame(h2o.predict(RF_Bal_200t,test_set.h2o))


RF_NoBal_ROC <- roc( MakeNumeric(RF_NoBal_Pred$predict), MakeNumeric(test_set$label) )
RF_Bal_ROC <- roc( MakeNumeric(RF_Bal_Pred$predict), MakeNumeric(test_set$label) )
RF_Bal_200t_ROC <- roc( MakeNumeric(RF_Bal_200t_Pred$predict), MakeNumeric(test_set$label) )



RF_NoBal_pred.auc <- prediction(MakeNumeric(RF_NoBal_Pred$predict), MakeNumeric(test_set$label))
RF_Bal_pred.auc <- prediction(MakeNumeric(RF_Bal_Pred$predict), MakeNumeric(test_set$label))
RF_Bal_200t_pred.auc <- prediction(MakeNumeric(RF_Bal_200t_Pred$predict), MakeNumeric(test_set$label))


ROCR::performance(RF_NoBal_pred.auc,'auc')
performance(RF_Bal_pred.auc,'auc')
performance(RF_Bal_200t_pred.auc,'auc')

pdf('plots/RFresults.pdf')
draw_confusion_matrix(confusionMatrix(data=RF_NoBal_Pred$predict, reference = test_set$label),'RF_NoBal', 'No','Yes',NULL,'rosybrown1')
draw_confusion_matrix(confusionMatrix(data=RF_Bal_Pred$predict, reference = test_set$label),'RF-Bal', 'No','Yes',NULL,'rosybrown1')
draw_confusion_matrix(confusionMatrix(data=RF_Bal_200t_Pred$predict, reference = test_set$label),'RF_Bal_200t', 'No','Yes',NULL,'rosybrown1')
dev.off()











