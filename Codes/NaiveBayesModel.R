Initialize()
# h2o.init(ip = 'localhost', port = 54321, nthreads= detectCores()-4, max_mem_size = '48g')

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

#### Naive-Bayes

nb1 <- naiveBayes(formula= label~., data=training_set)
NB1_pred <- predict(nb1, test_set)
confusionMatrix(data = NB1_pred , reference = test_set$label)
draw_confusion_matrix(confusionMatrix(data = NB1_pred , reference = test_set$label))

### h2o implementation:
features = colnames(training_set)[-ncol(training_set)]
label = 'label'
training_set.h2o = training_set
test_set.h2o = test_set
training_set.h2o$weight  <- ifelse(training_set$label=='YES',10,1)
test_set.h2o$weight <- ifelse(test_set$label=='YES', 10, 1)
training_set.h2o <- as.h2o(training_set)
test_set.h2o <- as.h2o(test_set)
nb1 = h2o.naiveBayes(x=features, y = label, training_frame = training_set.h2o, validation_frame = test_set.h2o, balance_classes = T, seed = 1398)
######



