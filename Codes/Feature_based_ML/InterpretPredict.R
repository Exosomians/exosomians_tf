source('Codes/Functions.R')
Initialize()
TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = F)
listOfSvmPreds <- readRDS('models/listOfSvmPreds.rds')



trainSet_lengthNormalized <-  readRDS('Data/trainSet_lengthNormalized.rds')
testSet_lengthNormalized <- readRDS('Data/testSet_lengthNormalized.rds')


listOfTitles <- c('SVM-weighted(1:10)-linearKernel', 
                  'SVM-weighted(1:10)-RadialKernel', 
                  'SVM-unWeighted-linearKernel', 
                  'SVM-weighted(1:10)-noStructFeature-linearKernel', 
                  'SVM-weighted(1:10)-noStructFeature-radialKernel',
                  'SVM-weighted(1:10)-SecondStruct(NO-Kmer)-linearKernel',
                  'SVM-weighted(1:10)-Minimal-linearKernel',
                  'SVM-weighted(1:10)-Min-SecStruct-linearKernel')

svmModelNames <- gsub('\\(1:10\\)|Kernel|ial|ear|uct|Feature|mal|ond','',
                      gsub('weighted|Weighted', 'Bal',listOfTitles))

### SVM predictions

SVMpreds <- lapply(listOfSvmPreds, as.character)
SVMpreds <- data.frame(do.call(cbind, SVMpreds))
colnames(SVMpreds) <- svmModelNames


### RF predictions

RFpreds <- lapply(listOfRFPreds, function(x) as.character(x$predict))
RFpreds <- data.frame(do.call(cbind, RFpreds))
colnames(RFpreds) <- names(listOfRFModels)


labelAtt <- TotalMatrixWithStruct[!TrainSetIndex,c('ic','ev', 'label','X0','X1','chr','length' )]






############# train label based feature distributions 

#TotalMatrixWithStruct <- subset(TotalMatrixWithStruct, select=-c(X))
trainSet_lenNorm <- trainSet_lengthNormalized[[4]]


Train_YES <- subset(trainSet_lenNorm, label=='YES' )
Train_NO <- subset(trainSet_lenNorm, label=='NO' )

NumericFeatureIndices <- sapply(1:ncol(trainSet_lenNorm), function(i) class(trainSet_lenNorm[,i]) %in% c('integer', 'numeric'))
NumericFeatures <- colnames(trainSet_lenNorm)[NumericFeatureIndices]

MakeTrainFeatureDf <- function(featureName){
  df = rbind(data.frame(feature=Train_YES[,featureName], name ='train.YES'), 
             data.frame(feature=Train_NO[,featureName], name = 'train.NO'))
  return(df)}


Train.evalFeature <- sapply(1:length(NumericFeatures), 
                            function(i) MakeTrainFeatureDf(NumericFeatures[i]), simplify = F )



pdf('plots/CheckFeatureDiscriminationQuality_lengthNormalized.pdf', height = 5)
TrainFeatureEvalPlots <- sapply(1:length(Train.evalFeature), function(i){
  
  p1=ggplot(Train.evalFeature[[i]], aes(color=name,x=feature))+
    geom_density()+theme_bw()+scale_fill_brewer()+
    labs(
      title = paste0(NumericFeatures[i], "-distribution"),
      subtitle = "Train Data label:NO, Train Data label:YES")+
    theme(
      plot.title = element_text(color = "dark blue", size = 10),
      plot.subtitle = element_text(color = "blue", size = 7)) + ylab('feature')
  
  
  p2=ggplot(Train.evalFeature[[i]], aes(x=name,y=feature))+
    geom_boxplot(aes(fill=name))+theme_bw()+scale_fill_brewer()+
    labs(
      title = paste0(NumericFeatures[i], "-distribution"),
      subtitle = "Train Data label:NO, Train Data label:YES")+
    theme(
      plot.title = element_text(color = "dark blue", size = 10),
      plot.subtitle = element_text(color = "blue", size = 7)) + ylab('feature')
  
  return(grid.arrange(p1, p2, nrow=1,ncol=2))}, simplify = F) 

dev.off()






############ Interpreting SVM predictions

svmBalRadPrected <- data.frame(svmPred = SVMpreds[,1], labelAtt)
trueExported_SVMwrongPredict <- subset(svmBalRadPrected, label=='YES' & svmPred=='NO')


MakeSVMFeatureDf <- function(featureName){
  df = rbind(data.frame(feature=Train_YES[,featureName], name ='train.YES'), 
             data.frame(feature=Train_NO[,featureName] ,name= 'train.NO'), 
             data.frame(feature=trueExported_SVMwrongPredict[,featureName], name='radBalSVM_wrongNO'))
  return(df)}


svm.evaLen <- MakeSVMFeatureDf('length')
svm.evalEV <- MakeSVMFeatureDf('ev')
svm.evalIC <- MakeSVMFeatureDf('ic')
svm.evalX0 <- MakeSVMFeatureDf('X0')
svm.evalX1 <- MakeSVMFeatureDf('X1')


p.svm_paired = ggplot(svm.evalX1, aes(x=name,y=feature+0.1))+geom_boxplot(aes(fill=name))+theme_bw()+scale_fill_brewer()
p.svm_paired = p.svm_paired + labs(
  title = "Paired-distribution",
  subtitle = "True label:YES, svmBalRad label:NO")+
  theme(
    plot.title = element_text(color = "dark blue", size = 10),
    plot.subtitle = element_text(color = "blue", size = 7),
    legend.position  = 'none') + scale_y_log10() + ylab('feature')



# pdf('plots/interpretBestSVM.pdf', height = 10)
# grid.arrange(p.svm_len, p.svm_ev, p.svm_ic, p.svm_unpaired, p.svm_paired, nrow=3, ncol=2)
# dev.off()

