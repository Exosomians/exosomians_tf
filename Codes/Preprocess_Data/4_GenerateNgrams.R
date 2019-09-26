#  In this file, we will try to find the ngrams(K-mers) 
#     of RNA sequences to use them as features. 
#     We'll make all the possible 4-grams and 2-grams
#     then, we'll check the count of the ngrams in each RNA sequence


# Required Libraries: parallel, reshape2, ggplot2
# library(parallel)
# library(reshape2)
# library(ggplot2)

options(stringsAsFactors = F)

DESIGN_MATRIX_PREP3 = 'Data/oldDataRefined/SecondStruct/3_DesignMat_SS_Label.csv'
DESIGN_MATRIX_PREP4 = 'Data/oldDataRefined/DesignMatrices/4_DesignMat_SS_Kmer_Label.csv'
NGRAM_MATRIX = 'Data/oldDataRefined/Ngrams/ngramMatrix.csv'

K_MER_SIZE = 4

source('Codes/Functions.R')
# Initialize()

designMatrix = read.csv(DESIGN_MATRIX_PREP3, stringsAsFactors = F)

All_Possible_Kmers <- Make_All_Possible_Kmers(K_MER_SIZE, 4, c('A','U','C','G'))

## Parallelize it using mclapply
# RNA_NucleotideKmers <- sapply(as.character(designMatrix$seq), makeKmerForNucleotide , K_MER_SIZE)
RNA_NucleotideKmers = mclapply(as.character(designMatrix$seq), makeKmerForNucleotide , K_MER_SIZE, mc.cores = detectCores()-2)
names(RNA_NucleotideKmers) = designMatrix$id

# summary(sapply(RNA_NucleotideKmers, length))

Nucleotide_KmerMatrix <- MakeFeatureSpecificMatrix(All_Possible_Kmers, RNA_NucleotideKmers, designMatrix$id)
# DrawFeatureDistribution(Nucleotide_KmerMatrix)
# summary(rowSums(Nucleotide_KmerMatrix))

write.csv(Nucleotide_KmerMatrix, NGRAM_MATRIX, quote = F)

designMatrixNgram <- cbind(designMatrix, Nucleotide_KmerMatrix)
write.csv(designMatrixNgram, DESIGN_MATRIX_PREP4, quote = F,row.names = F)