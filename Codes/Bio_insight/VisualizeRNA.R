if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("R4RNA")
library(R4RNA)


message("TRANSAT prediction in helix format")
transat_file <- system.file("extdata", "helix.txt", package = "R4RNA")
transat <- readHelix(transat_file)
message("RFAM structure in dot bracket format")
known_file <- system.file("extdata", "vienna.txt", package = "R4RNA")
known <- readVienna(known_file)
message("Work with basepairs instead of helices for more flexibility")
message("Breaks all helices into helices of length 1")

transat <- expandHelix(transat)
known <- expandHelix(known)


# simple arc plots
plotHelix(known, line = TRUE, arrow = TRUE)
mtext("Known Structure", side = 3, line = -2, adj = 0.1)



#Two structures for the same sequence can be visualized simultaneously, 
#allowing one to compare and contrast the two structures.
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)
mtext("TRANSAT\nPredicted\nStructure", side = 3, line = -5, adj = 0)
mtext("Known Structure", side = 1, line = -2, adj = 0)

message("Filter out helices above a certain p-value")
transat <- transat[which(transat$value <= 1e-3), ]


message("Assign colour to basepairs according to p-value")
transat$col <- col <- colourByValue(transat, log = TRUE)
message("Coloured encoded in 'col' column of transat structure")
plotDoubleHelix(transat, known, line = TRUE, arrow = TRUE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
          inset = 0.05, bty = "n", border = NA, cex = 0.75, title = "TRANSAT P-values")

plotOverlapHelix(transat, known, line = TRUE, arrow = TRUE, scale = FALSE)


message("Multiple sequence alignment of interest")
library(Biostrings)
fasta_file <- system.file("extdata", "fasta.txt", package = "R4RNA")
fasta <- as.character(readBStringSet(fasta_file))
message("Plot covariance in alignment")
plotCovariance(fasta, known, cex = 0.5)

plotCovariance(fasta, transat, cex = 0.5, conflict.col = "grey")

col <- colourByCovariation(known, fasta, get = TRUE)
plotCovariance(fasta, col, grid = TRUE, legend = FALSE)
legend("topright", legend = attr(col, "legend"), fill = attr(col, "fill"),
         inset = 0.1, bty = "n", border = NA, cex = 0.37, title = "Covariation")

col <- colourByUnknottedGroups(known, c("red", "blue"), get = TRUE)
plotCovariance(fasta, col, base.colour = TRUE, legend = FALSE, species = 23, grid = TRUE, text = TRUE)

