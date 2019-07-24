## In this file we write the sequences of our regions  
##  of interest in a fasta file, which is needed by 
##  ViennaRNA Package for RNA Secondary structure prediction

#   input: merged label and design matrix
#   output: fasta file


source('Codes/Functions.R')
Initialize()

TotalMatrix <- read.csv('Data/oldDataRefined/2_PrimaryDesignMat_label.csv',stringsAsFactors=F )
Sequence <- data.frame(id=TotalMatrix$id,seq=TotalMatrix$seq)
Sequence$seq <- as.character(Sequence$seq)


writeFasta<-function(data, filename){
  fastaLines = c()
  
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

head(Sequence)
writeFasta(Sequence, "Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta")



#################
# Secondary structure prediction using ViennaRNA tool

### bash 
# RNAfold  -i RegionsLenFilter.fasta > RNAseconStructPredict.txt
# grep '(\|>' RNAseconStructPredict.txt > RNAseconStructPredict_withoutSeq.txt

RNAseqFasta <- readRNAStringSet('Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta')
SeqID = names(RNAseqFasta)
RNAseqLength <- nchar(paste(RNAseqFasta))
hist(RNAseqLength,breaks = seq(1,520,10))

SecondaryStructureWithEnergy = readBStringSet("Data/oldDataRefined/SecondStruct/RNAseconStructPredict_withoutSeq.txt")
sum(names(SecondaryStructureWithEnergy)!=SeqID)

SecondaryStructureWithEnergy <- paste(SecondaryStructureWithEnergy)  
DotBracketSecondStruct  <-  substr(SecondaryStructureWithEnergy,1, RNAseqLength)

FreeEnergy <- gsub('[() ]','', 
                   substr(SecondaryStructureWithEnergy, RNAseqLength,nchar(SecondaryStructureWithEnergy)))

extraDotIndex <- substr(FreeEnergy,1,1)=='.'
FreeEnergy[extraDotIndex] = as.numeric(substr(FreeEnergy[extraDotIndex],2,
                                              nchar(FreeEnergy[extraDotIndex])))
FreeEnergy <- as.numeric(FreeEnergy)
class(FreeEnergy)





#################
## Generate Single-mer feature for dot-bracket 
All_Possible_DotBracket_SingleMer <- c('(','.')
DotBracket_SingleMer <- sapply(1:length(DotBracketSecondStruct),
                           function(i) unlist(str_split(gsub(')','(',DotBracketSecondStruct[i]),'')),simplify = F)

DotBracket_SingleMerMatrix <- MakeFeatureSpecificMatrix(All_Possible_DotBracket_SingleMer, DotBracket_SingleMer, Sequence$id)
head(DotBracket_SingleMerMatrix)
DrawFeatureDistribution(DotBracket_SingleMerMatrix)



################# IGNORE THIS SECTION -> UNINFORMATIVE FEATURE
### Generate Tri-mer feature for dot-bracket 

# add middle nucleotide to the K-mer 
All_Possible_DotBracket_3Kmer <-unlist(sapply(1:4, function(i) 
                                                paste0(c('A','U','C','G')[i],  Make_All_Possible_Kmers(3, 2,  c('(','.') )),simplify=F ))
# replace ')' with '('
DotBracket_3mer <- sapply(1:length(DotBracketSecondStruct), function(i) 
                            MakeKmerForDotBracket(Sequence$seq[i],
                                                  gsub( ')','(',DotBracketSecondStruct[i]), 3) )
# make design matrix for Kmer
DotBracket_3merMatrix <- MakeFeatureSpecificMatrix(All_Possible_DotBracket_3Kmer, DotBracket_3mer, Sequence$id)
DrawFeatureDistribution(DotBracket_3merMatrix)
#################


#### merging all design matrices together
RNAsecondaryStructFeatures <- cbind(FreeEnergy, DotBracket_SingleMerMatrix)
RNAsecondaryStructFeatures$DB <- DotBracketSecondStruct


RNAstructColNames <- colnames(RNAsecondaryStructFeatures)
RNAstructColNames <- sapply(RNAstructColNames, RefineSecondStructColNames, simplify = F)
colnames(RNAsecondaryStructFeatures) <- RNAstructColNames
RNAsecondaryStructFeatures$id <- rownames(RNAsecondaryStructFeatures)
rownames(RNAsecondaryStructFeatures) <- NULL

write.csv(RNAsecondaryStructFeatures, 'Data/oldDataRefined/SecondStruct/SecondStructFeatures_LenFilter.csv',row.names = F)
RNAsecondaryStructFeatures <- read.csv('Data/oldDataRefined/SecondStruct/SecondStructFeatures_LenFilter.csv', stringsAsFactors = F)





### merging the k-mer design matrix(filtered based on min-length=18) with secondary structures 
TotalMatrixWithStruct <- merge(TotalMatrix, RNAsecondaryStructFeatures, by.x='id',by.y='id',all.x=T)
write.csv(TotalMatrixWithStruct, 'Data/oldDataRefined/SecondStruct/3_DesignMat_SS_Label.csv',row.names = F)





## Additional notes for future
###################
# RNA binding proteins:

# check motif_scan > https://github.com/miha-skalic/motif_scan
# python2 /home/delaram/Tutorial/downloads/motif_scan/motif_scan.py -s AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCT > hits.tab
# deepbind

##################
# RNA motif annotation

# check biostars python code for motif extraction 
# https://github.com/cschu/biolib/blob/master/mdg_dt.py


##  convert dot_bracket notation > Ct
# bash
# RNAfold < ~/delaram/data/ExosomeProj/Data/RNAseconStructPredict.txt | b2ct > ct_file.ct

CtSecondStruct <- sapply(1:nrow(Sequence), 
                         function(i) makeCt(DotBracketSecondStruct[i], as.character(Sequence[i,'seq'])),
                         simplify = F)

names(CtSecondStruct) <- Sequence$id



###### Visualization
coord=ct2coord(ct)
RNAPlot(coord,hl=c("GGGUUU","AAAUUU"),seqcols=c(2,4),labTF=TRUE)
### Not working
bulge_loop(ct)
hairpin_loop(ct)
internal_loop(ct)
multi_branch_loop(ct)
stem(ct)
###################


