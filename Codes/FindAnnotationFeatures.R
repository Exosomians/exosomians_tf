## in this script, we'll try to find some secondary structure features
## such as number of hairpins, stemloops, ..., longest hairpin length, longest stem length

## for more information about these annotations, please check: https://viennarna.github.io/forgi/graph_tutorial.html


## import data 
data <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_element_string.csv', stringsAsFactors = F)
annotations <- data$element_string


motifCounts <- data.frame(do.call(rbind, lapply(annotations, .getNumberOfMotifs)))
colnames(motifCounts) <- c('fiveprime', 'threeprime', 'stem', 'interior_loop', 'multiloop', 'hairpin')

longestHairpins <- unlist(lapply(annotations, .getLongestMotif, 'h'))
longestStemLoop <- unlist(lapply(annotations, .getLongestMotif, 's'))

annotation_features <- data.frame(cbind(motifCounts, longestHairpins, longestStemLoop))
write.csv(annotation_features, 'Data/annotation_features.csv', row.names = F)


write.csv(cbind(data, annotation_features),'Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_PlusAnnot.csv', row.names=F)




.getNumberOfMotifs <- function(String){
  values <- replicate(6, 0)
  names(values) <-  c('f', 't', 's', 'i' , 'm', 'h' )
  for (i in 1:nchar(String)){
    if( substr(String, i, i) != substr(String, i-1, i-1) | i==1) values[substr(String,i,i)] = values[substr(String,i,i)]+1 
  } 
  return(values)
}




.getLongestMotif <- function(String, motif){
  i = 0 
  max_hairpin = 0 
  while (i < nchar(String)){
    counter = 0
    if (substr(String,i,i) == motif){
      counter = counter + 1
      for (j in i+1:nchar(String) ){
        if (substr(String,j,j) == motif ) counter = counter + 1
        else break
      }
      #print(paste0('inner ', motif ,' length is: ', counter))
      if(max_hairpin<counter) max_hairpin = counter
      i = i + counter
    }
    else i = i + 1
  }
  return(max_hairpin)
}

