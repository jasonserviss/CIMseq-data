packages <- c(
  "sp.scRNAseq",
  "sp.scRNAseqTesting",
  "tidyverse"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#load counts
#remove columns with "SI1.Singlet.1" in colnames. Early miSeq. Duplicated samples.
#several samples are only NA. find with apply(counts, 2, function(x) all(is.na(x)))
#Regev data identified by "SRR" in colnames
#singlets grepl("SRR", colnames(counts)) | grepl("Singlet", colnames(counts)) | grepl("NJA00102", colnames(counts)) | grepl("NJA00103", colnames(counts)) | grepl("NJA00104", colnames(counts)) | grepl("NJA00109", colnames(counts)) | grepl("NJA00204", colnames(counts)) | grepl("NJA00205", colnames(counts)) | grepl("NJA00206", colnames(counts))

path <- './inst/rawData/GFPmouse/count_table_180112.txt'
counts <- read.table(path, header = TRUE, sep = "\t")
counts <- counts[, !grepl("SI1.Singlet.1", colnames(counts))]

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}
#move genes to rownames
rownames(data) <- data$HGN
data$HGN <- NULL

#label singlets and multiplets
#martin already named these so I am just switching to the naming convention
#I typically use.
mul <- stringr::str_detect(colnames(data), "Doublet")
colnames(data) <- ifelse(
  mul,
  paste("m.", colnames(data), sep = ""),
  paste("s.", colnames(data), sep = "")
)
colnames(data) <- stringr::str_replace(colnames(data), "\\.Singlet", "")
colnames(data) <- stringr::str_replace(colnames(data), "\\.Doublet", "")

#remove "htseq" suffix
colnames(data) <- stringr::str_replace(colnames(data), "(.*)\\.htseq$", "\\1")

#extract ERCC
ercc <- str_detect(rownames(data), "^ERCC\\-[0-9]*$")
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
