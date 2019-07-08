source('Codes/Functions.R')
Initialize()

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_PlusAnnot.csv', stringsAsFactors = F)
NumericFeatures <- c("fiveprime",  "threeprime", "stem", "interior_loop",
                     "multiloop", "hairpin", "longestHairpins", "longestStemLoop")


Train_YES <- subset(TotalMatrixWithStruct, label=='YES')
Train_NO <- subset(TotalMatrixWithStruct, label=='NO')

Train_YES[, NumericFeatures] <- Train_YES[, NumericFeatures]/Train_YES$length
Train_NO[, NumericFeatures] <- Train_NO[, NumericFeatures]/Train_NO$length



MakeTrainFeatureDf <- function(featureName){
  df = rbind(data.frame(feature=Train_YES[,featureName], name ='train.YES'), 
             data.frame(feature=Train_NO[,featureName], name = 'train.NO'))
  return(df)}


Train.evalFeature <- sapply(1:length(NumericFeatures), 
                            function(i) MakeTrainFeatureDf(NumericFeatures[i]), simplify = F )




pdf('plots/CheckSecondStructFeature_lengthNormalized.pdf', height = 5)
TrainFeatureEvalPlots <- sapply(1:length(Train.evalFeature), function(i){
  
  p1=ggplot(Train.evalFeature[[i]], aes(color=name,x=feature))+
    geom_density()+theme_bw()+scale_fill_brewer()+
    labs(
      title = paste0(NumericFeatures[i], "-distribution"),
      subtitle = "length normalized")+
    theme(
      plot.title = element_text(color = "dark blue", size = 10),
      plot.subtitle = element_text(color = "blue", size = 7)) + ylab('feature')
  
  
  p2=ggplot(Train.evalFeature[[i]], aes(x=name,y=feature))+
    geom_boxplot(aes(fill=name))+theme_bw()+scale_fill_brewer()+
    labs(
      title = paste0(NumericFeatures[i], "-distribution"),
      subtitle = "length normalized")+
    theme(
      plot.title = element_text(color = "dark blue", size = 10),
      plot.subtitle = element_text(color = "blue", size = 7)) + ylab('feature')
  
  return(grid.arrange(p1, p2, nrow=1,ncol=2))}, simplify = F) 

dev.off()


