#run from package root
#source('inst/rawData/countsMgfp/counts_Mgfp.R')

packages <- c("sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

cat('Processing countsMgfpTumor.\n')
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download raw data
files <- c(
  'countsMgfpTumor_180810.txt',
  'countsMgfpTumorMeta_180810.txt'
)

paths <- file.path('./inst/rawData/countsMgfpTumor', files)

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

#label singlets and multiplets
#ids should include SINGLET plates only
ids <- c(
  "NJC01201\\.[A-Z][1-4]"
)
counts <- labelSingletsAndMultiplets(counts, ids)

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts, mincount = 0), ]

#remove low quality cells
lqc <- detectLowQualityCells(
  counts,
  geneName = "Actb",
  mincount = 1e4,
  quantileCut = 0.01
)
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#prepare metadata
plateData <- paths[grepl("Meta", paths)] %>%
  map(read_tsv) %>%
  bind_rows() %>%
  dplyr::mutate(sample = removeHTSEQsuffix(sample)) %>%
  dplyr::mutate(sample = labelSingletsAndMultiplets(sample, ids)) %>%
  annotatePlate(.) %>%
  annotateRow(.) %>%
  annotateColumn(.) %>%
  annotateMouse(
    .,
    plate = c("NJC01201"),
    mouse = c(12)
  ) %>%
  annotateTissue(
    .,
    plate = c("NJC01201"),
    tissue = c("colon")
  ) %>%
  dplyr::mutate(cellNumber = dplyr::if_else(
    stringr::str_detect(sample, "^s"), "Singlet", "Multiplet")
  ) %>%
  dplyr::mutate(filtered = dplyr::if_else(sample %in% colnames(counts), FALSE, TRUE))

#rename and save
countsMgfpTumor <- counts
countsMgfpTumorERCC <- countsERCC
countsMgfpTumorMeta <- plateData

if(all(!colnames(countsMgfpTumor) %in% countsMgfpTumorMeta$sample)) {
  stop("all counts data not present in meta data")
}
if(all(!colnames(countsMgfpTumorERCC) %in% countsMgfpTumorMeta$sample)) {
  stop("all ercc data not present in meta data")
}

save(
  countsMgfpTumor,
  countsMgfpTumorERCC,
  countsMgfpTumorMeta,
  file = "./data/countsMgfpTumor.rda",
  compress = "bzip2"
)
