## Installing H2O
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))

library(h2o)

h2o.init(ip = 'localhost', port = 54321, nthreads= detectCores()-4, max_mem_size = '48g')

h2o.clusterInfo()

exData = finalMat
colsToBeFactorized = c('id', 'label', 'chr', 'strand', 'annotation', 'rnaType')
exData[colsToBeFactorized] = lapply(exData[colsToBeFactorized], factor)

exData = as.h2o(exData)
str(exData)

?h2o.randomForest


trainSet = read.csv('Data/SecondaryDataFolded/Fold02_Train.csv')
testSet = read.csv('Data/SecondaryDataFolded/Fold02_Test.csv')

trainSet$weight = 0
trainSet$weight[trainSet$label=='Yes'] = 10
trainSet$weight[trainSet$label=='No'] = 1

testSet$weight = 0
testSet$weight[testSet$label=='Yes'] = 10
testSet$weight[testSet$label=='No'] = 1

colsToBeFactorized = c('id', 'label', 'chr', 'strand', 'annotation', 'rnaType')

trainSet[colsToBeFactorized] = lapply(trainSet[colsToBeFactorized], factor)
testSet[colsToBeFactorized] = lapply(testSet[colsToBeFactorized], factor)

trainSet = as.h2o(trainSet)
testSet = as.h2o(testSet)

features = c('chr', 'length', 'a', 'c', 'g', 't', 'strand', 'annotation', 'rnaType')
label = 'label'

rf = h2o.randomForest(x = features, y = label, training_frame = trainSet, model_id = 'Exosox1', validation_frame = testSet, nfolds = 5, balance_classes = T, ntrees = 50, max_depth = 20, seed = 1398)
rf2 = h2o.randomForest(x = features, y = label, training_frame = trainSet, model_id = 'Exosox1', validation_frame = testSet, nfolds = 5, balance_classes = T, ntrees = 200, max_depth = 20, seed = 1398)
rf3 = h2o.randomForest(x = c(features, 'seq'), y = label, training_frame = trainSet, model_id = 'Exosox1', validation_frame = testSet, nfolds = 5, balance_classes = T, ntrees = 50, max_depth = 20, seed = 1398)
rf4 = h2o.randomForest(x = features, y = label, training_frame = trainSet, model_id = 'Exosox1', validation_frame = testSet, nfolds = 5, balance_classes = F, ntrees = 50, max_depth = 20, seed = 1398, weights_column = 'weight')
