#run from package root
#source('./inst/rawData/countsSorted2/countsSorted2_171116.R')

library(sp.scRNAseqData)
library(dplyr)
cat('Processing countsSorted2.\n')

googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download data
googledrive::drive_download(
  file = 'countsSorted2_171116.txt',
  path = './inst/rawData/countsSorted2/countsSorted2_171116.txt',
  overwrite = TRUE
)

googledrive::drive_download(
  file = 'countsSortedMeta2_171116.txt',
  path = './inst/rawData/countsSorted2/countsSortedMeta2_171116.txt',
  overwrite = TRUE
)

#load counts
path <- './inst/rawData/countsSorted2/countsSorted2_171116.txt'
counts <- read.table(path, header = TRUE)

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
lqc.totalCounts <- detectLowQualityCells.totalCounts(counts, mincount = 4e4)
lqc.housekeeping <- detectLowQualityCells.housekeeping(counts, geneName = "Actb", quantileCut = 0.01)
lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(counts, countsERCC, percentile = 0.99)
lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac
print(paste0("Removing a total of ", sum(!lqc), " cells based on the calculated metrics."))
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#setup metadata
plateData <- loadMetaData('./inst/rawData/countsSorted2/countsSortedMeta2_171116.txt') %>%
dplyr::mutate(sample = removeHTSEQsuffix(sample)) %>%
dplyr::mutate(sample = labelSingletsAndMultiplets(sample, "NJB00201")) %>%
annotatePlate(.) %>%
annotateRow(.) %>%
annotateColumn(.) %>%
annotateCellNumber(., plate = "NJB00201", row = 1:8, column = 1:12) %>%
countsSorted2_cellTypes(.) %>%
dplyr::mutate(filtered = dplyr::if_else(sample %in% colnames(counts), FALSE, TRUE))

#rename and save
if(!all(colnames(counts) %in% plateData$sample)) {
  stop("all counts data not present in meta data")
}
if(!all(colnames(countsERCC) %in% plateData$sample)) {
  stop("all ercc data not present in meta data")
}

countsSorted2 <- counts
countsSortedERCC2 <- countsERCC
countsSortedMeta2 <- plateData
save(
  countsSorted2,
  countsSortedERCC2,
  countsSortedMeta2,
  file = "./data/countsSorted2.rda",
  compress = "bzip2"
)
