## In this file we write the sequences of our regions  
##  of interest in a fasta file, which is needed by 
##  ViennaRNA Package for RNA Secondary structure prediction

#   input: merged label and design matrix
#   output: fasta file

# Required Libraries: parallel, Biostrings, seqinr, stringr, reshape2, ggplot2
# library(parallel)
# library(Biostrings)
# library(seqinr)
# library(stringr)
# library(reshape2)
# library(ggplot2)

options(stringsAsFactors = F)

source('Codes/Functions.R')
# Initialize()

SEED = 1
THREADS = detectCores()-2

DESIGN_MATRIX_PREP2 = 'Data/oldDataRefined/DesignMatrices/2_all_PrimaryDesignMat_label_db.csv'
DESIGN_MATRIX_PREP3='Data/oldDataRefined/SecondStruct/3_DesignMat_SS_Label.csv'
SEQUENCES_FASTA = 'Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta'
SECONDARY_STRUCTURES = 'Data/oldDataRefined/SecondStruct/RNAseconStructPredict.txt'
SECONDARY_STRUCTURES_WO_SEQ= 'Data/oldDataRefined/SecondStruct/RNAseconStructPredict_withoutSeq.txt'
SECONDARY_STRUCTURES_LEN_FILTER='Data/oldDataRefined/SecondStruct/SecondStructFeatures_LenFilter.csv'

TCGA_SEQ_MIN_LENGTH = 20
TCGA_SEQ_MAX_LENGTH = 200

TotalMatrix <- read.csv(DESIGN_MATRIX_PREP2)

# Sampling the data
# All H, All S, Sample T in a way that #S=#T
numOfSerumSamples = length(which(TotalMatrix$db == 'S'))
tcgaSamples = subset(TotalMatrix, db=='T')


summary(tcgaSamples$length)
# TCGA data are miRNA-seq, so their sequences length are small.

# Filtering TCGA sequences by length,
# 20 < length < 200 is desired.
tcgaFilteredByLength = subset(tcgaSamples, length>= TCGA_SEQ_MIN_LENGTH &
                                length<=TCGA_SEQ_MAX_LENGTH)

# Removing sequences that have N.
tcgaWhichSequencesContainsN = mclapply(strsplit(tcgaFilteredByLength$seq, ''),
                                       function(eachSequence) 'N'%in%eachSequence,
                                       mc.cores = detectCores()-2)
tcgaWhichSequencesContainsN = which(unlist(tcgaWhichSequencesContainsN))

tcgaFilteredDueToHavingN = tcgaFilteredByLength[-tcgaWhichSequencesContainsN, ]

# Data is imbalanced (Label NO is much more than YES),
# We downsample the TCGA sequences and extract 105374 sequences from them,
# 105374 is the number of Serum (exRNA) sequences which we assumme them as YES labels.
set.seed(SEED)
tcgaDownsampledIndexes = sample(seq(nrow(tcgaFilteredDueToHavingN)), numOfSerumSamples)
tcgaDownsampled = tcgaFilteredByLength[tcgaDownsampledIndexes, ]

tcgaSamples = tcgaDownsampled
attr(tcgaSamples, 'description') =
  'TCGA Samples, Filtered by length, Filtered due to having N in the sequence, Downsampled to the number of Serum samples'
# attributes(tcgaSamples)


rm(tcgaFilteredByLength,
   tcgaWhichSequencesContainsN,
   tcgaFilteredDueToHavingN,
   tcgaDownsampledIndexes,
   tcgaDownsampled)
gc()

designMatrixBalanced = rbind(tcgaSamples, subset(TotalMatrix, db!='T'))

# Make the "id" of sequences more informative,
# Add the "db" (H, S, T) and the "label" (Y, N) to the "id".
designMatrixBalanced$id = paste0(designMatrixBalanced$id,
                                 ifelse(designMatrixBalanced$label=='YES', 'Y', 'N'),
                                 designMatrixBalanced$db)

designMatrix = designMatrixBalanced
attr(designMatrix, 'description') =
  'Design Matrix of the sequences from Hani, Serum , and TCGA databases; balanced based on label.'

rm(TotalMatrix,
   tcgaSamples,
   designMatrixBalanced)
gc()

### Refining the sequences:
# Sequences on the + strand should have been complemented,
# but all sequences are have been reverseComplemented,
# so the sequences on the + strand should be reveresed.
positiveStrandSequences = designMatrix[grepl('[+]', designMatrix$id), ]$seq
class(positiveStrandSequences)
positiveStrandSequences = RNAStringSet(positiveStrandSequences)
positiveStrandSequences = reverse(positiveStrandSequences)
positiveStrandSequences = as.character(positiveStrandSequences)
designMatrix[grepl('[+]', designMatrix$id), ]$seq = positiveStrandSequences

rm(positiveStrandSequences)

seqinr::write.fasta(sequences = strsplit(designMatrix$seq, ''),
                    names = designMatrix$id,
                    file.out = SEQUENCES_FASTA)


#################
# Secondary structure prediction + Free Energy using ViennaRNA tool

### Installing Vienna-RNA on Ubuntu 18.04
# wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb
# sudo apt install libgsl23 libgslcblas0
# sudo dpkg -i viennarna_2.4.14-1_amd64.deb
# sudo apt install -f

### Bash
# JOBS=`nproc`
# SEQUENCES_FASTA=Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta
# SECONDARY_STRUCTURES=Data/oldDataRefined/SecondStruct/RNAseconStructPredict.txt
# SECONDARY_STRUCTURES_WO_SEQ=Data/oldDataRefined/SecondStruct/RNAseconStructPredict_withoutSeq.txt
# RNAfold --noPS --jobs=$JOBS -i $SEQUENCES_FASTA > $SECONDARY_STRUCTURES
# grep '(\|>' $SECONDARY_STRUCTURES > $SECONDARY_STRUCTURES_WO_SEQ

runRNAfoldCommand =
  sprintf("JOBS=%s;
          SEQUENCES_FASTA=%s;
          SECONDARY_STRUCTURES=%s;
          SECONDARY_STRUCTURES_WO_SEQ=%s;
          RNAfold  --noPS --jobs=$JOBS -i $SEQUENCES_FASTA > $SECONDARY_STRUCTURES;
          grep '(\\|>' $SECONDARY_STRUCTURES > $SECONDARY_STRUCTURES_WO_SEQ",
          THREADS,
          SEQUENCES_FASTA,
          SECONDARY_STRUCTURES,
          SECONDARY_STRUCTURES_WO_SEQ)
runRNAfoldCommand = gsub('\n', '', runRNAfoldCommand)


system(runRNAfoldCommand)


## Parsing ViennaRNA Output:
## Extracting "Secondary Structure" of "Dot Bracket" form,
## and "Free Energy" values.

# (ViennaRNA predicts the secondary structures of the RNA sequences;
# In addition, ViennaRNA calculates/estiamtes the free energy of each sequence.)

RNAseqFasta = readRNAStringSet(SEQUENCES_FASTA)
SeqID = names(RNAseqFasta)
RNAseqLength <- nchar(paste(RNAseqFasta))
# hist(RNAseqLength,breaks = seq(1,520,5), xlim = c(1,100))

SecondaryStructureWithEnergy = readBStringSet(SECONDARY_STRUCTURES_WO_SEQ)
sum(names(!(SecondaryStructureWithEnergy)%in%SeqID))

SecondaryStructureWithEnergy <- paste(SecondaryStructureWithEnergy)  
DotBracketSecondStruct  <-  substr(SecondaryStructureWithEnergy,1, RNAseqLength)

FreeEnergy <- gsub('[() ]','', 
                   substr(SecondaryStructureWithEnergy, RNAseqLength,nchar(SecondaryStructureWithEnergy)))

extraDotIndex <- substr(FreeEnergy,1,1)=='.'
FreeEnergy[extraDotIndex] = as.numeric(substr(FreeEnergy[extraDotIndex],2,
                                              nchar(FreeEnergy[extraDotIndex])))
FreeEnergy <- as.numeric(FreeEnergy)




## Generate 1-mer features ['.', '('] from Dot-Bracket sequences,
# Dot-Bracket form of secondary structures were generated by ViennaRNA
All_Possible_DotBracket_SingleMer <- c('(','.')
DotBracket_SingleMer <- sapply(1:length(DotBracketSecondStruct),
                               function(i) unlist(str_split(gsub(')','(',DotBracketSecondStruct[i]),'')),simplify = F)

DotBracket_SingleMerMatrix <- MakeFeatureSpecificMatrix(All_Possible_DotBracket_SingleMer, DotBracket_SingleMer, designMatrix$id)
# head(DotBracket_SingleMerMatrix)
# DrawFeatureDistribution(DotBracket_SingleMerMatrix)


### IGNORE THIS SECTION -> UNINFORMATIVE FEATURE
## Generate Tri-mer feature for dot-bracket 
# # add middle nucleotide to the K-mer 
# All_Possible_DotBracket_3Kmer <-unlist(sapply(1:4, function(i) 
#   paste0(c('A','U','C','G')[i],  Make_All_Possible_Kmers(3, 2,  c('(','.') )),simplify=F ))
# # replace ')' with '('
# DotBracket_3mer <- sapply(1:length(DotBracketSecondStruct), function(i) 
#   MakeKmerForDotBracket(Sequence$seq[i],
#                         gsub( ')','(',DotBracketSecondStruct[i]), 3) )
# # make design matrix for Kmer
# DotBracket_3merMatrix <- MakeFeatureSpecificMatrix(All_Possible_DotBracket_3Kmer, DotBracket_3mer, Sequence$id)
# DrawFeatureDistribution(DotBracket_3merMatrix)
###

## Merging all design matrices together
RNAsecondaryStructFeatures <- cbind(FreeEnergy, DotBracket_SingleMerMatrix)
RNAsecondaryStructFeatures$DB <- DotBracketSecondStruct

RNAstructColNames <- colnames(RNAsecondaryStructFeatures)
RNAstructColNames <- sapply(RNAstructColNames, RefineSecondStructColNames, simplify = F)
colnames(RNAsecondaryStructFeatures) <- RNAstructColNames
RNAsecondaryStructFeatures$id <- rownames(RNAsecondaryStructFeatures)
rownames(RNAsecondaryStructFeatures) <- NULL

write.csv(RNAsecondaryStructFeatures, SECONDARY_STRUCTURES_LEN_FILTER, row.names = F)
# RNAsecondaryStructFeatures <- read.csv(SECONDARY_STRUCTURES_LEN_FILTER, stringsAsFactors = F)


## Merging Generated Features
### merging the design matrix (filtered based on min-length=18) with secondary structures
TotalMatrixWithStruct <- merge(designMatrix, RNAsecondaryStructFeatures, by.x='id',by.y='id',all.x=T)
write.csv(TotalMatrixWithStruct, DESIGN_MATRIX_PREP3, row.names = F)


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
## CT : Connectivity Table
# Bash on Ubuntu 18.04
# SECONDARY_STRUCTURES=Data/oldDataRefined/SecondStruct/RNAseconStructPredict.txt
# SECONDARY_STRUCTURES_CT=Data/oldDataRefined/SecondStruct/ct_file.ct
# RNAfold < $SECONDARY_STRUCTURES | b2ct > $SECONDARY_STRUCTURES_CT


# CtSecondStruct <- sapply(1:nrow(Sequence), 
#                          function(i) makeCt(DotBracketSecondStruct[i], as.character(Sequence[i,'seq'])),
#                          simplify = F)
# 
# names(CtSecondStruct) <- Sequence$id
# 
# 
# 
# ###### Visualization
# coord=ct2coord(ct)
# RNAPlot(coord,hl=c("GGGUUU","AAAUUU"),seqcols=c(2,4),labTF=TRUE)
# ### Not working
# bulge_loop(ct)
# hairpin_loop(ct)
# internal_loop(ct)
# multi_branch_loop(ct)
# stem(ct)
# ###################