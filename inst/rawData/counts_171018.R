#NJB00101 is the singlet plate, NJB00103 is doublets (according to your scheme).

path <- './inst/rawData/counts_171018.txt'
countsMe <- read.table(path, sep = "\t", header = TRUE)

#move genes to rownames
rownames(countsMe) <- countsMe$HGN
countsMe$HGN <- NULL

#annotate singlets and multiplets
s <- grepl("NJB00101", colnames(countsMe))
colnames(countsMe) <- ifelse(
    s,
    paste("s.", colnames(countsMe), sep = ""),
    paste("m.", colnames(countsMe), sep = "")
)

#remove "htseq" suffix
colnames(countsMe) <- gsub("(.*)\\.htseq$", "\\1", colnames(countsMe))

#extract ERCC
ercc <- grepl("^ERCC\\-[0-9]*$", rownames(countsMe))
countsERCC <- countsMe[ercc, ]

if(dim(countsERCC)[1] != 92) {
    stop("Couldn't detect all ERCC reads.")
}

countsMe <- countsMe[!ercc, ]

#remove "bad" genes
countsMe <- countsMe[rowSums(countsMe) > 0, ]

nonGenes <- c(
"__no_feature", "__ambiguous", "__too_low_aQual",
"__not_aligned", "__alignment_not_unique"
)

countsMe <- countsMe[!rownames(countsMe) %in% nonGenes, ]

#remove cells with poor coverage
mincount <- 1e5
countsERCC <- countsERCC[, colSums(countsMe) > mincount]
countsMe <- countsMe[, colSums(countsMe) > mincount]

countsMe <- as.matrix(countsMe)
countsERCC <- as.matrix(countsERCC)

#rename and save
countsSorted1 <- countsMe
countsSortedERCC1 <- countsERCC
save(
    countsSorted1,
    countsSortedERCC1,
    file = "./data/countsSorted1.rda",
    compress = "bzip2"
)

#source('./inst/rawData/counts_171018.R')