---
title: "Exosomians"
output: html_notebook
---

```{r}
## R
options(stringsAsFactors = F)

source('Codes/Functions.R')
Initialize()

SEED = 1
THREADS = detectCores()-2
```


# Preprocessing

## Data Importion
**Importing raw reads (FASTQ files)**

## Read Alignment
**Aligning reads to reference assembly (hg38)**

Code: bash
Tools: hisat2

## Region Extraction
Code: bash
Tools: bedtools

## Region Labeling
Code: bash
Tools: bedtools

## Sequence Extraction and Refinement
Code: R
Packages: BSgenome.Hsapiens.UCSC.hg38
```{r}
# R
# Required Libraries: parallel, Biostrings

## Sampling the data
## All H, All S, Sample T in a way that #S=#T
numOfSerumSamples = length(which(designMatrix$db == 'S'))
tcgaSamples = subset(designMatrix, db=='T')

## TCGA data are miRNA-seq, so their sequences length are small.
summary(tcgaSamples$length)

## Filtering TCGA sequences by length,
## 20 < length < 200 is desired.
tcgaFilteredByLength = subset(tcgaSamples, length>= TCGA_SEQ_MIN_LENGTH &
                                length<=TCGA_SEQ_MAX_LENGTH)

## Removing sequences that have N.
tcgaWhichSequencesContainsN = mclapply(strsplit(tcgaFilteredByLength$seq, ''),
                                       function(eachSequence) 'N'%in%eachSequence,
                                       mc.cores = detectCores()-2)
tcgaWhichSequencesContainsN = which(unlist(tcgaWhichSequencesContainsN))

tcgaFilteredDueToHavingN = tcgaFilteredByLength[-tcgaWhichSequencesContainsN, ]
```

## Label Balancing

```{r}
## R
TCGA_SEQ_MIN_LENGTH = 20
TCGA_SEQ_MAX_LENGTH = 200

## Data is imbalanced (Label NO is much more than YES),
## We downsample the TCGA sequences and extract 105374 sequences from them,
## 105374 is the number of Serum (exRNA) sequences which we assumme them as YES labels.
set.seed(SEED)
tcgaDownsampledIndexes = sample(seq(nrow(tcgaFilteredDueToHavingN)), numOfSerumSamples)
tcgaDownsampled = tcgaFilteredByLength[tcgaDownsampledIndexes, ]

tcgaSamples = tcgaDownsampled
attr(tcgaSamples, 'description') =
  'TCGA Samples, Filtered by length, Filtered due to having N in the sequence, Downsampled to the number of Serum samples'
# attributes(tcgaSamples)


designMatrixBalanced = rbind(tcgaSamples, subset(TotalMatrix, db!='T'))

## Make the "id" of sequences more informative,
## Add the "db" (H, S, T) and the "label" (Y, N) to the "id".
designMatrixBalanced$id = paste0(designMatrixBalanced$id,
                                 ifelse(designMatrixBalanced$label=='YES', 'Y', 'N'),
                                 designMatrixBalanced$db)

designMatrix = designMatrixBalanced
attr(designMatrix, 'description') =
  'Design Matrix of the sequences from Hani, Serum , and TCGA databases; balanced based on label.'


rm(tcgaFilteredByLength,
   tcgaWhichSequencesContainsN,
   tcgaFilteredDueToHavingN,
   tcgaDownsampledIndexes,
   tcgaDownsampled,
   tcgaSamples,
   designMatrixBalanced)
gc()
```

# Processing (Feature Generation)

## General
### Chromosome

**Extract "Chromosome Number" of each sequence as "chr"**

### Cytoband

**Extract "Cytoband Location" of each sequence as "cyto"**

### Length

**Extract "Length" of each sequence as "len"**

## Primary Sequence

### ACGU Sequence

**Extract "ACGU Sequence" of each sequence as "seq"**

(This feature has been extracted in "Preprocessing: Sequence Extraction" phase)

### 1-mer (A/C/G/U)

**Extract "Number of A/C/G/U" in each sequence as "a/c/g/u"**

### 4-mer (AAAA/AAAC/AAAG/...)

**Extract "Number of Each 4-mer" in each sequence as "[4-MER_SEQUENCE]"**

```{r}
## R
## Required Libraries: parallel

K_MER_SIZE = 4

# options(stringsAsFactors = F)
# source('Codes/Functions.R')
# Initialize()

designMatrix

All_Possible_Kmers <- Make_All_Possible_Kmers(K_MER_SIZE, 4, c('A','U','C','G'))

## Parallelize it using mclapply
# RNA_NucleotideKmers <- sapply(as.character(designMatrix$seq), makeKmerForNucleotide , K_MER_SIZE)
RNA_NucleotideKmers = mclapply(as.character(designMatrix$seq), makeKmerForNucleotide , K_MER_SIZE,
                               mc.cores = detectCores()-2)
names(RNA_NucleotideKmers) = designMatrix$id

Nucleotide_KmerMatrix <- MakeFeatureSpecificMatrix(All_Possible_Kmers, RNA_NucleotideKmers, designMatrix$id)

designMatrix <- cbind(designMatrix, Nucleotide_KmerMatrix)
```


### DeepBind Scores

**Extract "DeepBind Scores" of each sequence as "[DEEPBIND_MODEL_NAME]"**

```{bash}
## Bash
## Downloading and installing DeepBind
cd ~/R/Projects/exosomians; mkdir Apps; cd Apps
wget http://tools.genes.toronto.edu/deepbind/deepbind-v0.11-linux.tgz
tar xf deepbind-v0.11-linux.tgz
```

```{r}
## R
## Extracting all DeepBind models name related to Homo Sapiens.
## ./deepbind/db/db.tsv

deepbindModels = read.delim('Apps/deepbind/db/db.tsv', comment.char = '#')
# unique(deepbindModels$Species) # "Homo sapiens"

deepbindModels = deepbindModels[deepbindModels$Species=='Homo sapiens', ]
deepbindModelsId = deepbindModels$ID

deepbindParams = list.files('Apps/deepbind/db/params/')
deepbindParams = gsub('.txt', '', deepbindParams)

deepbindAvailableModels = deepbindModelsId[deepbindModelsId%in%deepbindParams]

## Writing DeepBind model ids
write.table(deepbindAvailableModels, file = 'Apps/deepbind/models.ids', 
            quote = F, row.names = F, col.names = F)

## Writing sequences
write.table(designMatrix$seq, file = 'Apps/deepbind/sequences.seq', 
            quote = F, row.names = F, col.names = F)
```

```{bash}
## Bash
## Calculating DeepBind scores

## Single-threaded: Not feasible for many sequences
# system('cd Apps/deepbind; ./deepbind models.ids < sequences.seq > deepbind.output')

## Mutli-threaded: Ghavi!,
## Parallelize it by splitting models id into several parts,
## and run deepbind on multicores

cd ~/R/Projects/exosomians/Apps/deepbind
mkdir models; cp models.ids models/
THREADS=`nproc`; THREADS=`expr $THREADS - 2`
split --additional-suffix=".ids" --number=l/$THREADS models.ids
rm models.ids

## (CAUTION: Run the code below in a TMUX environment! It might take a while...)
for model in *.ids; do ../deepbind $model < ../sequences.seq > ${model/.ids/.output} & done
paste *.output | column -s $'\t' -tn > all.output
```

```{r}
## R
deepbindScores = read.table('Apps/deepbind/models/all.output', header = T)
designMatrix <- cbind(designMatrix, deepbindScores)

deepbindDictionary = deepbindModels[deepbindModels$ID %in% deepbindAvailableModels, ]
write.csv(deepbindDictionary, 'Apps/deepbind/deepbind.dict', quote = F, row.names = F)
```

## Secondary Sequence

### Dot-bracket Sequence

**Extract "Dot-bracket Sequence" of each sequence as "DB"**

```{r}
## R
## Required Libraries: parallel, Biostrings, seqinr

SEQUENCES_FASTA = 'Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta'

write.fasta(sequences = strsplit(designMatrix$seq, ''),
            names = designMatrix$id,
            file.out = SEQUENCES_FASTA)
```

```{bash}
## Bash
## Secondary structure prediction + Free Energy using ViennaRNA tool
## Installing Vienna-RNA on Ubuntu 18.04

wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb
sudo apt install libgsl23 libgslcblas0
sudo dpkg -i viennarna_2.4.14-1_amd64.deb
sudo apt install -f

JOBS=`nproc`
SEQUENCES_FASTA=Data/oldDataRefined/SecondStruct/RegionsLenFilter.fasta
SECONDARY_STRUCTURES=Data/oldDataRefined/SecondStruct/RNAseconStructPredict.txt
SECONDARY_STRUCTURES_WO_SEQ=Data/oldDataRefined/SecondStruct/RNAseconStructPredict_withoutSeq.txt
RNAfold  --noPS --jobs=$JOBS -i $SEQUENCES_FASTA > $SECONDARY_STRUCTURES;
grep '(\\|>' $SECONDARY_STRUCTURES > $SECONDARY_STRUCTURES_WO_SEQ",
```

```{r}
## R
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
```

### Free Energy

**Extract "Free Energy" of each sequence as "FreeEnergy"**

```{r}
## R
FreeEnergy <- gsub('[() ]','', 
                   substr(SecondaryStructureWithEnergy, RNAseqLength,nchar(SecondaryStructureWithEnergy)))

extraDotIndex <- substr(FreeEnergy,1,1)=='.'
FreeEnergy[extraDotIndex] = as.numeric(substr(FreeEnergy[extraDotIndex],2,
                                              nchar(FreeEnergy[extraDotIndex])))
FreeEnergy <- as.numeric(FreeEnergy)
```

### 1-mer

**Extract "Number of ./(" in each sequence as "X0/X1"**

```{r}
## R
All_Possible_DotBracket_SingleMer <- c('(','.')
DotBracket_SingleMer <- sapply(1:length(DotBracketSecondStruct),
                               function(i) unlist(str_split(gsub(')','(',DotBracketSecondStruct[i]),'')),
                               simplify = F)

DotBracket_SingleMerMatrix <- MakeFeatureSpecificMatrix(All_Possible_DotBracket_SingleMer, DotBracket_SingleMer, designMatrix$id)
# DrawFeatureDistribution(DotBracket_SingleMerMatrix)
```
```{r}
## Merging all design matrices together
RNAsecondaryStructFeatures <- cbind(FreeEnergy, DotBracket_SingleMerMatrix)
RNAsecondaryStructFeatures$DB <- DotBracketSecondStruct

RNAstructColNames <- colnames(RNAsecondaryStructFeatures)
RNAstructColNames <- sapply(RNAstructColNames, RefineSecondStructColNames, simplify = F)
colnames(RNAsecondaryStructFeatures) <- RNAstructColNames
RNAsecondaryStructFeatures$id <- rownames(RNAsecondaryStructFeatures)
rownames(RNAsecondaryStructFeatures) <- NULL

## Merging the design matrix (filtered based on min-length=18) with secondary structures
TotalMatrixWithStruct <- merge(designMatrix, RNAsecondaryStructFeatures, by.x='id',by.y='id',all.x=T)
write.csv(TotalMatrixWithStruct, DESIGN_MATRIX_PREP3, row.names = F)
```