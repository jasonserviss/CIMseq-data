#run from package root
#source('inst/rawData/countsSorted1/counts_171018.R')
#NJB00101 is the singlets plate, NJB00103 is doublets (according to your scheme).

path <- './inst/rawData/sortedMultiplets_171018/counts_171018.txt'
countsMe <- read.table(path, sep = "\t", header = TRUE)

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
lqc <- detectLowQualityCells(countsMe)
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
