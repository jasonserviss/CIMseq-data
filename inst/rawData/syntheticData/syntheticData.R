##THIS IS NOT CURRENTLY FINISHED
library(sp.scRNAseqTesting)

outPath <- './data'

tmp <- syntheticMultuplets(
  nGenes = 2000,
  nCells = 100,
  nCellTypes = 10,
  perplexity = 10,
  target = 20,
  singletExpansion = 20
)

singlets <- tmp[[1]]
multuplets <- tmp[[2]]
uObj <- tmp[[3]]

colnames(multuplets) <- gsub(
  "^([A-Z0-9]*)\\..*",
  "\\1",
  colnames(multuplets)
)

table <- calculateConnections(testMultuplets, type = "multuplets")

names <- c(
    paste("s", colnames(singlets), sep = "."),
    paste("m", colnames(multuplets), sep = ".")
)

syntheticData <- cbind(singlets, multuplets)
colnames(syntheticData) <- names
syntheticData <- as.matrix(syntheticData)

#rename and save
syntheticDataUnsupervised <- uObj
syntheticDataTable <- table

save(
  syntheticData,
  syntheticDataUnsupervised,
  syntheticDataTable,
  file = paste(outPath, "syntheticData.rda", sep = "/"),
  compress = "bzip2"
)

#source('./inst/rawData/syntheticData.R')

