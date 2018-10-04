#run from package root
#source('./inst/rawData/countsSorted2/countsSorted2_171116.R')

packages <- c("sp.scRNAseqData", "EngeMetadata")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

cat('Processing countsSorted2.\n')

googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")
basePath <- 'data/Sorted.multiplets'

#download data
googledrive::drive_download(
  file = file.path(basePath, 'raw_data/countsSorted2_171116.txt'),
  path = './inst/rawData/countsSorted2/countsSorted2_171116.txt',
  overwrite = TRUE
)

plates <- c("NJB00204", "NJB00201")
plateData <- purrr::map_dfr(plates, function(p) {
  metadata(p, file.path(basePath, 'annotation'))
}) %>%
  dplyr::mutate(prefix = dplyr::if_else(cellNumber == "Singlet", "s.", "m.")) %>%
  dplyr::mutate(sample = paste0(prefix, unique_key, ".", Well)) %>%
  dplyr::select(-prefix) %>%
  dplyr::select(sample, dplyr::everything())

#load counts
paths <- './inst/rawData/countsSorted2/countsSorted2_171116.txt'
counts <- read.table(paths, header = TRUE)

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
lqc.housekeeping <- detectLowQualityCells.housekeeping(counts, geneName = "ACTB", quantileCut = 0.01)
lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(counts, countsERCC, percentile = 0.99)
lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac
print(paste0("Removing a total of ", sum(!lqc), " cells based on the calculated metrics."))
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#setup metadata
plateData <- plateData %>%
  dplyr::mutate(filtered = dplyr::if_else(sample %in% colnames(counts), FALSE, TRUE))


if(!all(colnames(counts) %in% plateData$sample)) {
  stop("all counts data not present in meta data")
}
if(!all(colnames(countsERCC) %in% plateData$sample)) {
  stop("all ercc data not present in meta data")
}

#rename and save as .rda
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

#save processed data as text and upload to drive
metaPath <- './inst/rawData/countsSorted2/countsSortedMeta2.txt'
countsPath <- './inst/rawData/countsSorted2/countsSorted2.txt'
erccPath <- './inst/rawData/countsSorted2/countsSortedERCC2.txt'

write_tsv(countsMgfpMeta, path = metaPath)
write_tsv(as.data.frame(countsMgfp), path = countsPath)
write_tsv(as.data.frame(countsMgfpERCC), path = erccPath)

googledrive::drive_upload(metaPath, file.path(basePath, "processed_data/countsSortedMeta2.txt"))
googledrive::drive_upload(countsPath, file.path(basePath, "processed_data/countsSorted2.txt"))
googledrive::drive_upload(erccPath, file.path(basePath, "processed_data/countsSortedERCC2.txt"))

#delete all text files
trash <- map(c(paths, metaPath, countsPath, erccPath), file.remove)
