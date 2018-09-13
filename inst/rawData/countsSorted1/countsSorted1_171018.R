#run from package root
#source('inst/rawData/countsSorted1/countsSorted1_171018.R')

cat('Processing countsSorted1.\n')
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download data
googledrive::drive_download(file = 'countsSorted1_171018.txt', path = './inst/rawData/countsSorted1/countsSorted1_171018.txt', overwrite = TRUE)

#NJB00101 is the singlets plate, NJB00103 is doublets (according to your scheme).

path <- './inst/rawData/countsSorted1/countsSorted1_171018.txt'
counts <- read.table(path, sep = "\t", header = TRUE)

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#annotate singlets and multiplets
counts <- labelSingletsAndMultiplets(counts, "NJB00101")

#remove "htseq" suffix
counts <- removeHTSEQsuffix(counts)

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove genes with low counts
counts <- counts[detectLowQualityGenes(counts), ]

#remove cells with poor coverage
lqc.totalCounts <- detectLowQualityCells.totalCounts(counts, mincount = 4e4)
lqc.housekeeping <- detectLowQualityCells.housekeeping(counts, geneName = "Actb", quantileCut = 0.01)
lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(counts, countsERCC, percentile = 0.99)
lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac
print(paste0("Removing a total of ", sum(!lqc), " cells based on the calculated metrics."))
countsERCC <- countsERCC[, lqc]
counts <- counts[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#rename and save
countsSorted1 <- counts
countsSortedERCC1 <- countsERCC
save(
  countsSorted1,
  countsSortedERCC1,
  file = "./data/countsSorted1.rda",
  compress = "bzip2"
)
