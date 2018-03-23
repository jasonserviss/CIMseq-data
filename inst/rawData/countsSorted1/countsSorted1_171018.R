#run from package root
#source('inst/rawData/countsSorted1/counts_171018.R')

library(googledrive)
gs_auth(token = "data/googlesheets_token.rds")

#download data
drive_download(file = 'countsSorted1_171018.txt', path = './inst/rawData/countsSorted1/countsSorted1_171018.txt', overwrite = TRUE)

#NJB00101 is the singlets plate, NJB00103 is doublets (according to your scheme).

path <- './inst/rawData/countsSorted1/countsSorted1_171018.txt'
countsMe <- read.table(path, sep = "\t", header = TRUE)

#check for NAs
if(sum(is.na(countsMe)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
countsMe <- moveGenesToRownames(countsMe)

#annotate singlets and multiplets
countsMe <- labelSingletsAndMultiplets(countsMe, "NJB00101")

#remove "htseq" suffix
countsMe <- removeHTSEQsuffix(countsMe)

#extract ERCC
ercc <- detectERCCreads(countsMe)
countsERCC <- countsMe[ercc, ]
countsMe <- countsMe[!ercc, ]

#remove non-genes
countsMe <- countsMe[!detectNonGenes(countsMe), ]

#remove genes with low counts
countsMe <- countsMe[detectLowQualityGenes(countsMe), ]

#remove cells with poor coverage
lqc <- detectLowQualityCells(
  countsMe,
  mincount = 4e4,
  quantileCut = 0.01
)
countsERCC <- countsERCC[, lqc]
countsMe <- countsMe[, lqc]

#coerce to matrix
countsMe <- convertCountsToMatrix(countsMe)
countsERCC <- convertCountsToMatrix(countsERCC)

#rename and save
countsSorted1 <- countsMe
countsSortedERCC1 <- countsERCC
save(
    countsSorted1,
    countsSortedERCC1,
    file = "./data/countsSorted1.rda",
    compress = "bzip2"
)
