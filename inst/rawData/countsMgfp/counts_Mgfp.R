#run from package root
#source('inst/rawData/countsMgfp/counts_Mgfp.R')

packages <- c("sp.scRNAseqData", "tidyverse", "EngeMetadata")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

cat('Processing countsMgfp.\n')

basePath <- 'data/Mouse.gut.architecture'
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download metadata
plates <- c(
  "NJA00101", "NJA00102", "NJA00103", "NJA00104", "NJA00107", "NJA00109",
  "NJA00110", "NJA00111", "NJA00201", "NJA00204", "NJA00205", "NJA00206",
  "NJA00402", "NJA00403", "NJA00404", "NJA00405", "NJA00406", "NJA00408",
  "NJA00409", "NJA00411", "NJA00412", "NJA00413", "NJA00602", "NJA00608",
  "NJA00609", "NJA00801", "NJA01202", "NJA01203", "NJA01205", "NJA01301",
  "NJA01302", "NJA01303", "NJA01401"
)

plateData <- purrr::map_dfr(plates, function(p) {
  path <- file.path(
    basePath, 'annotation',
    stringr::str_replace(p, "(.{6}).*", "\\1")
  )
  metadata(p, path)
}) %>%
  dplyr::mutate(prefix = dplyr::if_else(cellNumber == "Singlet", "s.", "m.")) %>%
  dplyr::mutate(sample = paste0(prefix, unique_key, ".", Well)) %>%
  dplyr::select(-prefix) %>%
  dplyr::select(sample, dplyr::everything())

#download raw counts data
files <- c(
  'countsMgfp_180316.txt',
  'countsMgfp_180409.txt',
  'countsMgfp_180427.txt',
  'countsMgfp_180604.txt',
  'countsMgfp_180810.txt',
  'countsMgfp_180910.txt'
)

paths <- file.path('./inst/rawData/countsMgfp', files)
trash <- map2(files, paths, function(file, path) {
  googledrive::drive_download(file, path, overwrite = TRUE)
})
rm(trash)

#load counts
dataFiles <- paths[!grepl("Meta", paths)]
loaded <- map(dataFiles, read.table, header = TRUE, sep = "\t")
counts <- reduce(loaded, full_join, by = "HGN")

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#remame samples with old nomenclature
colnames(counts) <- renameMgfpSamples(colnames(counts))

#label singlets and multiplets
#ids should include SINGLET plates only
singlets <- plateData[plateData$cellNumber == "Singlet", "sample"][[1]]
singlets <- gsub("^..(.*)", "\\1", singlets)
counts <- labelSingletsAndMultiplets(counts, singlets)
namesPreFilter <- colnames(counts)

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

#update metadata to show filtered samples
plateData <- mutate(plateData, filtered = if_else(sample %in% colnames(counts), FALSE, TRUE))

#rename and save
countsMgfp <- counts
countsMgfpERCC <- countsERCC
countsMgfpMeta <- plateData

#check all count samples in meta and vice versa
c1 <- all(!countsMgfpMeta$sample %in% namesPreFilter)
c2 <- all(!colnames(countsMgfp) %in% countsMgfpMeta$sample)
if(c1 & c2) {
  stop("all counts data not present in meta data")
}
if(all(!colnames(countsMgfpERCC) %in% countsMgfpMeta$sample)) {
  stop("all ercc data not present in meta data")
}

#save as .rda
save(
  countsMgfp,
  countsMgfpERCC,
  countsMgfpMeta,
  file = "./data/countsMgfp.rda",
  compress = "bzip2"
)

#save processed data as text and upload to drive
metaPath <- './inst/rawData/countsMgfp/countsMgfpMeta.txt'
countsPath <- './inst/rawData/countsMgfp/countsMgfp.txt'
erccPath <- './inst/rawData/countsMgfp/countsMgfpERCC.txt'

write_tsv(countsMgfpMeta, path = metaPath)
write_tsv(as.data.frame(countsMgfp), path = countsPath)
write_tsv(as.data.frame(countsMgfpERCC), path = erccPath)

googledrive::drive_upload(metaPath, file.path(basePath, "processed_data/countsMgfpMeta.txt"))
googledrive::drive_upload(countsPath, file.path(basePath, "processed_data/countsMgfp.txt"))
googledrive::drive_upload(erccPath, file.path(basePath, "processed_data/countsMgfpERCC.txt"))

#delete all text files
trash <- map(c(paths, metaPath, countsPath, erccPath), file.remove)
