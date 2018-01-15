packages <- c(
    "sp.scRNAseq",
    "sp.scRNAseqTesting",
    "tidyverse"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)


#load counts (171116 are prefixed with NJB00201 and NJB00204)
path <- './inst/rawData/counts_171116.txt'
data <- read.table(path, header = TRUE)

#move genes to rownames
rownames(data) <- data$HGN
data$HGN <- NULL

#subset relevant samples corresponding to date
data <- data[, str_detect(colnames(data), "NJB00201") | str_detect(colnames(data), "NJB00204")]

#label singlets and multiplets
sng <- str_detect(colnames(data), "NJB00201")
colnames(data) <- ifelse(
    sng,
    paste("s.", colnames(data), sep = ""),
    paste("m.", colnames(data), sep = "")
)

#remove "htseq" suffix
colnames(data) <- str_replace(colnames(data), "(.*)\\.htseq$", "\\1")

#extract ERCC
ercc <- grepl("^ERCC\\-[0-9]*$", rownames(data))
countsERCC <- data[ercc, ]

if(dim(countsERCC)[1] != 92) {
    stop("Couldn't detect all ERCC reads.")
}

data <- data[!ercc, ]

#remove "bad" genes
data <- data[rowSums(data) > 0, ]

nonGenes <- c(
  "__no_feature", "__ambiguous", "__too_low_aQual",
  "__not_aligned", "__alignment_not_unique"
)

data <- data[!rownames(data) %in% nonGenes, ]

#remove cells with poor coverage
mincount <- 1e5
countsERCC <- countsERCC[, colSums(data) > mincount]
data <- data[, colSums(data) > mincount]

data <- as.matrix(data)
countsERCC <- as.matrix(countsERCC)

#rename and save
countsSorted2 <- data
countsSortedERCC2 <- countsERCC
save(
  countsSorted2,
  countsSortedERCC2,
  file = "./data/countsSorted2.rda",
  compress = "bzip2"
)

#source('./inst/rawData/counts_171116.R')
