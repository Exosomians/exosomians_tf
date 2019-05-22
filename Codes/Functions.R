
Initialize = function()
{
  #### Loading required libraries ####
  options(stringsAsFactors = F)
  
  #BiocManager::install('pacman')
  
  pac <- list( 'ggplot2', 'limma', 'pheatmap','Biostrings','RNAstructureModuleMiner','RRNA',
              'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer', 'caret',
              'plyr', 'gridExtra','caTools', 'h2o','gtools','stringr')
  
  print(paste(pac , lapply(pac, require, character.only = TRUE), sep = ' : '))
  
  pacman::p_load( 'ggplot2', 'limma', 'pheatmap', 'Biostrings','RNAstructureModuleMiner','RRNA',
                 'VennDiagram','e1071' ,'reshape2', 'ggrepel', 'RColorBrewer', 'caret',
                 'plyr','gridExtra','caTools','h2o','gtools','stringr')
}



#### Functions ####
FindRepetitionIndex = function(aVector)
{
  i = 1
  repetitionIndex = integer(length(aVector))
  repetitionIndex[!duplicated(aVector)] = i
  aVector[!duplicated(aVector)] = NA
  
  while (!all(is.na(aVector))) {
    i = i+1
    repetitionIndex[!duplicated(aVector) & !is.na(aVector)] = i
    aVector[!duplicated(aVector)] = NA
  }
  
  repetitionIndex  
}



CalcSummaryForCellLines = function(selectedCounts)
{
  cellLinesName = unique(gsub(colnames(selectedCounts), pattern = '_.*', replacement =  ''))
  cellLinesSummary = lapply(cellLinesName, function(aCellLineName)
  {
    theCellLineCounts = selectedCounts[, grepl(colnames(selectedCounts), pattern = aCellLineName)]
    summary(rowMeans(theCellLineCounts))
  })
  names(cellLinesSummary) = cellLinesName
  
  cellLinesSummary
}



makeKmerForNucleotide <- function(String, KmerSize){
  sapply( 1:(nchar(String) - KmerSize+1 ), 
          function(i)
            substr(String, i, i+KmerSize-1 )
  )
}



MakeKmerForDotBracket <- function(nucleotide, dotbracket, kmerSize){
  sapply( 1:(nchar(dotbracket) - kmerSize+1 ), 
          function(i)
            paste0(substr(nucleotide,i+1,i+1),substr(dotbracket, i, i+kmerSize-1 ))
  )
}



MakeFeatureSpecificMatrix <- function(All_possible_features, extractedFromData, id){
  
  MatchData2Features <- sapply(1:length(extractedFromData), 
                               function(i) table(factor(extractedFromData[[i]], levels = All_possible_features )) ,simplify = F)
  
  FeatureMatrix <- as.data.frame(do.call(rbind, MatchData2Features))
  colnames(FeatureMatrix) <- All_possible_features
  rownames(FeatureMatrix) <- id
  return(FeatureMatrix)
}




Make_All_Possible_Kmers <- function(KmerSize, VectorSize ,VectorOfElements){
  Permutes <- gtools::permutations(n=VectorSize , r=KmerSize , v= VectorOfElements ,repeats=T)
  Permutes <- apply( Permutes, 1, function(x) paste0(x, collapse=''))
  return(Permutes)
}




DrawFeatureDistribution <- function(FeatureMat ){
  featureDis <- melt(colSums(FeatureMat))
  featureDis <- data.frame(kmer=rownames(featureDis),Count=featureDis$value)
  ggplot(featureDis,aes(x=kmer,y=Count,color='black'))+
    geom_bar(stat = 'identity',color="dark blue", fill="cadetblue2",width=0.6)+xlab('Feature')+theme_bw()+coord_flip()
}




DrawFeatureDistribution_LabelBased <- function(FeatureMat){
  FeatureMatSplit <-  split(FeatureMat, FeatureMat$label)
  FeatureMatSplit <- lapply(FeatureMatSplit, function(x) melt(colSums(x[,1:ncol(x)-1])))
  FeatureMatSplit <- lapply(FeatureMatSplit, function(x) x/sum(x))
  FeatureMatSplit <- sapply(1:length(FeatureMatSplit), function(i)
    data.frame(Feature=rownames(FeatureMatSplit[[i]]),
               Count=FeatureMatSplit[[i]]$value, 
               label=names(FeatureMatSplit)[i]),simplify = F)
  
  FeatureMat <- do.call(rbind, FeatureMatSplit)
  ggplot(FeatureMat,aes(x=Feature, y=Count, fill=label))+
    geom_bar(stat = 'identity',width=0.6, position = position_dodge(), color='black')+
    xlab('Feature')+ylab('Frequency')+theme_bw()#+coord_flip()
}






RefineSecondStructColNames <- function(StructColName){
  Characters <- unlist(str_split( StructColName ,'') )
  Characters[Characters=='.'] <- '0'
  Characters[Characters=='('] <- '1'
  return(paste(Characters,collapse = ''))
}



MakeNumeric <- function(ChrToNum){ifelse(as.character(ChrToNum)=='YES',1,0)}



SplitLabelFromFeatures <- function(aMatrix){
  Splited <- list(subset(aMatrix,select=-label),  aMatrix$label)
  names(Splited) <- c('features', 'label')
  return(Splited)
}


getLabel <- function( aSplitedList ){ aSplitedList[['label']] }
getFeatures <- function( aSplitedList ){ aSplitedList[['features']]}


svmModel <- function( training_set, Kernel, Weights){
  svm(formula= label~.,
      data = training_set , 
      type ='C-classification', 
      kernel = Kernel ,
      class.weights = Weights, 
      probability = T, 
      cachesize = 10000)
}



draw_confusion_matrix <- function(cm, ModelName ,
                                  FirstClassName=NULL, SecondClassName=NULL,
                                  FirstClassColor=NULL, SecondClassColor=NULL) {
  
  if(is.null(FirstClassName)){FirstClassName='Class1'}
  if(is.null(SecondClassName)){SecondClassName='Class1'}
  
  if(is.null(FirstClassColor)){FirstClassColor='#3F97D0'}
  if(is.null(SecondClassColor)){SecondClassColor='darkgoldenrod1'}

  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(paste0(ModelName,'- Confusion Matrix'), cex.main=1)
  
  # create the matrix 
  rect(150, 430, 240, 370, col=FirstClassColor) # class1 color 
  text(195, 435, FirstClassName , cex=1.2)
  rect(250, 430, 340, 370, col=SecondClassColor) # calss2 color
  text(295, 435, SecondClassName, cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=SecondClassColor)
  rect(250, 305, 340, 365, col=FirstClassColor)
  text(140, 400, FirstClassName , cex=1.2, srt=90)
  text(140, 335, SecondClassName, cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

