## in this script, we'll try to find some secondary structure features
## such as number of hairpins, stemloops, ..., longest hairpin length, longest stem length

## for more information about these annotations, please check: https://viennarna.github.io/forgi/graph_tutorial.html

# Required Libraries: 

options(stringsAsFactors = F)

SEED = 1
# THREADS = detectCores()-2

# DESIGN_MATRIX_PREP5 = 'Data/oldDataRefined/DesignMatrices/5_DesignMat_SS_Kmer_Deepbind_Label.csv'
DESIGN_MATRIX_PREP6 = 'Data/oldDataRefined/DesignMatrices/6_DesignMat_SS_Kmer_Deepbind_Forgi_Label.csv'
DESIGN_MATRIX_PREP7 = 'Data/oldDataRefined/DesignMatrices/7_DesignMat_SS_Kmer_Deepbind_ForgiAnnoted_Label.csv'

SECONDARY_STRUCTURES_BULGE_GRAPH = 'Data/oldDataRefined/SecondStruct/annotation_features.csv'

source('Codes/Functions.R')
Initialize()

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


## Import data 
designMatrix <- read.csv(DESIGN_MATRIX_PREP6, stringsAsFactors = F)
annotations <- designMatrix$element_string

motifCounts <- data.frame(do.call(rbind, lapply(annotations, .getNumberOfMotifs)))
colnames(motifCounts) <- c('fiveprime', 'threeprime', 'stem', 'interior_loop', 'multiloop', 'hairpin')

longestHairpins <- unlist(lapply(annotations, .getLongestMotif, 'h'))
longestStemLoop <- unlist(lapply(annotations, .getLongestMotif, 's'))

annotation_features <- data.frame(cbind(motifCounts, longestHairpins, longestStemLoop))
write.csv(annotation_features, SECONDARY_STRUCTURES_BULGE_GRAPH, row.names = F)

write.csv(cbind(subset(designMatrix,select=-element_string_number), annotation_features),
          DESIGN_MATRIX_PREP7, row.names=F)


