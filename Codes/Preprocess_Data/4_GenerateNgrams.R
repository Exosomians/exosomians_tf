#  In this file, we will try to find the ngrams(K-mers) 
#     of RNA sequences to use them as features. 
#     We'll make all the possible 4-grams and 2-grams
#     then, we'll check the count of the ngrams in each RNA sequence



source('Codes/Functions.R')
Initialize()
designMat = read.csv('Data/oldDataRefined/DesignMatrices/3_DesignMat_SS_Label.csv', stringsAsFactors = F)


K_MER_SIZE = 4
All_Possible_Kmers <- Make_All_Possible_Kmers(K_MER_SIZE, 4, c('A','U','C','G'))


RNA_NucleotideKmers <- sapply(as.character(designMat$seq), makeKmerForNucleotide , K_MER_SIZE)
summary(sapply(RNA_NucleotideKmers, length))


Nucleotide_KmerMatrix <- MakeFeatureSpecificMatrix(All_Possible_Kmers, RNA_NucleotideKmers, designMat$id)
DrawFeatureDistribution(Nucleotide_KmerMatrix)

summary(rowSums(Nucleotide_KmerMatrix))
write.csv(Nucleotide_KmerMatrix,'Data/oldDataRefined/Ngrams/ngramMatrix.csv',quote = F)

designMatNgram <- cbind(designMat, Nucleotide_KmerMatrix)
write.csv(designMatNgram,'Data/oldDataRefined/DesignMatrices/4_DesignMat_SS_Kmer_Label.csv',quote = F,row.names = F)


