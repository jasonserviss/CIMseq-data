#run from package root
#source('./inst/rawData/countsRegev/counts_Regev.R')
#GSE92332

library(sp.scRNAseqData)
cat('Processing countsRegev.\n')
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download data
googledrive::drive_download(
  file = 'counts_Regev.txt',
  path = './inst/rawData/countsRegev/counts_Regev.txt',
  overwrite = TRUE
)

#load counts
paths <- './inst/rawData/countsRegev/counts_Regev.txt'
counts <- read.table(paths, header = TRUE, sep = "\t")

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#remove "htseq" suffix
counts <- removeHTSEQsuffix(counts)

#All samples are singlets so label them as such.
counts <- labelSingletsAndMultiplets(counts, "SRR")

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts), ]

#remove low quality cells
lqc.totalCounts <- detectLowQualityCells.totalCounts(counts, mincount = 1e4)
lqc.housekeeping <- detectLowQualityCells.housekeeping(counts, geneName = "Actb", quantileCut = 0.01)
lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(counts, countsERCC, percentile = 0.99)
lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac
print(paste0("Removing a total of ", sum(!lqc), " cells based on the calculated metrics."))
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#rename and save
countsRegev <- counts
countsRegevERCC <- countsERCC

save(
  countsRegev,
  countsRegevERCC,
  file = "./data/countsRegev.rda",
  compress = "bzip2"
)

countsPath <- './inst/rawData/countsSorted2/countsRegev.txt'
erccPath <- './inst/rawData/countsSorted2/countsRegevERCC.txt'

write_tsv(as.data.frame(countsMgfp), path = countsPath)
write_tsv(as.data.frame(countsMgfpERCC), path = erccPath)

googledrive::drive_upload(countsPath, file.path(basePath, "processed_data/countsRegev.txt"))
googledrive::drive_upload(erccPath, file.path(basePath, "processed_data/countsRegevERCC.txt"))

#delete all text files
trash <- map(c(paths, countsPath, erccPath), file.remove)
