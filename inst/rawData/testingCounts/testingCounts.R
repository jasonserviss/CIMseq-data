#run from package root
#source('./inst/rawData/testingCounts/testingCounts.R')
#generates a data.frame of counts that resemble a small subset of real counts
#data to be used for testing.

#generate "normal" counts
nCellTypes <- 2
nCells <- 5
nGenes <- 10

synth <- sapply(1:nCellTypes, function(x) {
  set.seed(x + 8923023)
  rnbinom(nGenes * nCells, mu = 2^runif(nGenes, 0, 5), size = x)
})

counts <- matrix(as.numeric(synth), nrow = nGenes)

#add low quality cell
set.seed(989023)
counts <- cbind(
  counts,
  matrix(rnbinom(nGenes, mu = 2^runif(nGenes, 0, 1), size = 1), ncol = 1)
)

#add low quality gene
set.seed(83049)
counts <- rbind(
  counts,
  matrix(rnbinom((nGenes * 3) + 3, mu = 2^runif(nGenes, 0, 1), size = 1), nrow = 3)
)

#add names
colnames(counts) <- paste0(LETTERS[1:ncol(counts)], ".htseq")
rownames(counts) <- c(letters[1:11], "ERCC-1", "__alignment_not_unique")

#convert to data frame
counts <- as.data.frame(counts)

#Move gene names (rownames) to a column (instead of rownames) as is typical
#in HTSeq output
counts <- cbind(HGN = rownames(counts), counts)
rownames(counts) <- 1:nrow(counts)

#rename and save
testingCounts <- counts
save(testingCounts, file = "./data/testingCounts.rda", compress = "bzip2")

