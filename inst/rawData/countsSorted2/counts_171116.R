#run from package root
#source('./inst/rawData/countsSorted2/counts_171116.R')

#load counts (171116 are prefixed with NJB00201 and NJB00204)
path <- './inst/rawData/countsSorted2/counts_171116.txt'
counts <- read.table(path, header = TRUE)
bool1 <- grepl("NJB00201", colnames(counts)) | grepl("NJB00204", colnames(counts))
bool2 <- colnames(counts) == "HGN"
counts <- counts[, bool1 | bool2]

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#label singlets and multiplets
counts <- labelSingletsAndMultiplets(counts, "NJB00201")

#remove "htseq" suffix
counts <- removeHTSEQsuffix(counts)

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts), ]

#remove cells with poor coverage
lqc <- detectLowQualityCells(
  counts,
  mincount = 4e4,
  quantileCut = 0.01
)
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#rename and save
countsSorted2 <- counts
countsSortedERCC2 <- countsERCC
save(
  countsSorted2,
  countsSortedERCC2,
  file = "./data/countsSorted2.rda",
  compress = "bzip2"
)
