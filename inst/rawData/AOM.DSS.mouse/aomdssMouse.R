#run from package root
#source('inst/rawData/countsMgfp/counts_Mgfp.R')

packages <- c("sp.scRNAseqData", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

cat('Processing aomdssMouse\n')
basePath <- 'data/AOM.DSS.mouse'
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download metadata
plates <- c("NJC01201")

plateData <- purrr::map_dfr(plates, function(p) {
  path <- file.path(basePath, 'annotation')
  metadata(p, path)
}) %>%
  dplyr::mutate(prefix = dplyr::if_else(cellNumber == "Singlet", "s.", "m.")) %>%
  dplyr::mutate(sample = paste0(prefix, unique_key, ".", Well)) %>%
  dplyr::select(-prefix) %>%
  dplyr::select(sample, dplyr::everything())

#download raw data
files <- c('countsMgfpTumor_180810.txt')
paths <- file.path('./inst/rawData/AOM.DSS.mouse', files)

trash <- purrr::map2(files, paths, function(file, path) {
  googledrive::drive_download(file, path, overwrite = TRUE)
})
rm(trash)

#load counts
dataFiles <- paths[!grepl("Meta", paths)]
loaded <- purrr::map(dataFiles, read.table, header = TRUE, sep = "\t")
counts <- reduce(loaded, dplyr::full_join, by = "HGN")

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#label singlets and multiplets
#ids should include SINGLET plates only
singlets <- plateData[plateData$cellNumber == "Singlet", "sample"][[1]]
singlets <- gsub("^..(.*)", "\\1", singlets)
counts <- labelSingletsAndMultiplets(counts, singlets)

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts, mincount = 0), ]

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

#add filtered column to metadata
plateData <- plateData %>%
  dplyr::mutate(filtered = dplyr::if_else(sample %in% colnames(counts), FALSE, TRUE))

#rename and save
countsAomdss <- counts
countsAomdssERCC <- countsERCC
countsAomdssMeta <- plateData

if(all(!colnames(countsAomdss) %in% countsAomdssMeta$sample)) {
  stop("all counts data not present in meta data")
}
if(all(!colnames(countsAomdssERCC) %in% countsAomdssMeta$sample)) {
  stop("all ercc data not present in meta data")
}

#save as .rda
save(
  countsAomdss,
  countsAomdssERCC,
  countsAomdssMeta,
  file = "./data/countsAomdss.rda",
  compress = "bzip2"
)

#save processed data as text and upload to drive
metaPath <- './inst/rawData/countsMgfp/countsAomdssMeta.txt'
countsPath <- './inst/rawData/countsMgfp/countsAomdss.txt'
erccPath <- './inst/rawData/countsMgfp/countsAomdssERCC.txt'

write_tsv(countsAomdssMeta, path = metaPath)
write_tsv(as.data.frame(countsAomdss), path = countsPath)
write_tsv(as.data.frame(countsAomdssERCC), path = erccPath)

googledrive::drive_upload(metaPath, file.path(basePath, "processed_data/countsAomdssMeta.txt"))
googledrive::drive_upload(countsPath, file.path(basePath, "processed_data/countsAomdss.txt"))
googledrive::drive_upload(erccPath, file.path(basePath, "processed_data/countsAomdssERCC.txt"))

#delete all text files
trash <- map(c(paths, metaPath, countsPath, erccPath), file.remove)

