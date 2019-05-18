#  In this file, we will try to find the ngrams(K-mers) 
#     of RNA sequences to use them as features. 
#     We'll make all the possible 4-grams and 2-grams
#     then, we'll check the count of the ngrams in each RNA sequence


designMat = read.csv('Data/PrimaryDesignMat.csv')

source('Codes/Functions.R')
Initialize()

K_MER_SIZE = 4
All_Possible_Kmers <- Make_All_Possible_Kmers(K_MER_SIZE, 4, c('A','T','C','G'))



RNA_NucleotideKmers <- sapply(as.character(designMat$seq), makeKmerForNucleotide , K_MER_SIZE)
summary(sapply(RNA_NucleotideKmers, length))


Nucleotide_KmerMatrix <- MakeFeatureSpecificMatrix(All_Possible_Kmers, RNA_NucleotideKmers, designMat$id)
DrawFeatureDistribution(Nucleotide_KmerMatrix)

summary(rowSums(Nucleotide_KmerMatrix))
write.csv(Nucleotide_KmerMatrix,'Data/ngramMatrix.csv',quote = F)

designMatNgram <- cbind(designMat, Nucleotide_KmerMatrix)
write.csv(designMatNgram,'Data/NgramLen4DesignMatrix.csv',quote = F)












