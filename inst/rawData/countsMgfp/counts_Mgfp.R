#run from package root
#source('inst/rawData/countsMgfp/counts_Mgfp.R')

packages <- c("sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

cat('Processing countsMgfp.\n')
googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

#download raw data
files <- c(
  'countsMgfp_180316.txt',
  'countsMgfpMeta_180316.txt',
  'countsMgfp_180409.txt',
  'countsMgfpMeta_180409.txt',
  'countsMgfp_180427.txt',
  'countsMgfpMeta_180427.txt',
  'countsMgfp_180604.txt',
  'countsMgfpMeta_180604.txt',
  'countsMgfp_180810.txt',
  'countsMgfpMeta_180810.txt',
  'countsMgfp_180910.txt',
  'countsMgfpMeta_180910.txt'
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

#label singlets and multiplets
#ids should include SINGLET plates only
ids <- c(
  "Singlet", "NJA00102", "NJA00103", "NJA00104",
  "NJA00109", "NJA00204", "NJA00205", "NJA00206",
  "NJA00402", "NJA00403", "NJA00411", "NJA00412",
  "NJA00404", "NJA00405", "NJA00408", "NJA00409",
  "NJA00602", "NJA00608", "NJA00609", "NJA00801",
  "NJA01202\\.[A-Z][0-1][0-8]", "NJA01202\\.[A-Z]09",
  "NJA01203\\.[A-Z][0-1][4-9]", "NJA01203\\.[A-Z]1[0-1]",
  "NJA01205\\.[A-Z][1-4]", "NJA01301\\.[A-Z][0][1-9]",
  "NJA01301\\.[A-Z][1][0-8]", "NJA01302\\.[A-Z][0][1-9]",
  "NJA01302\\.[A-Z][1][0-8]", "NJA01303\\.[A-Z][0][1-9]",
  "NJA01303\\.[A-Z][1][0-8]", "NJA01401\\.[A-Z][0][1-9]",
  "NJA01401\\.[A-Z][1][0-8]"
)
counts <- labelSingletsAndMultiplets(counts, ids)

#remame samples with old nomenclature
colnames(counts) <- renameMgfpSamples(colnames(counts))

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

#prepare metadata
plateData <- paths[grepl("Meta", paths)] %>%
  map(read_tsv) %>%
  bind_rows() %>%
  dplyr::mutate(sample = removeHTSEQsuffix(sample)) %>%
  dplyr::mutate(sample = labelSingletsAndMultiplets(sample, ids))%>%
  dplyr::mutate(sample = renameMgfpSamples(sample)) %>%
  annotatePlate(.) %>%
  annotateRow(.) %>%
  annotateColumn(.) %>%
  annotateGFP(
    .,
    plate = list(
      c("NJA00111"), c("NJA00110"), c("NJA00103"),
      c("NJA00104"), c("NJA00201"), c("NJA00110")
    ),
    row = list(1:8, 1:8, 1:8, 1:8, 1:8, 1:8),
    column = list(1:12, 1:12, 1:6, 1:6, 1:12, 7:12)
  ) %>%
  dplyr::mutate(GFP = dplyr::if_else(
    plate %in% c(
      "NJA00204", "NJA00205", "NJA00206", "NJA00402", "NJA00403", "NJA00411",
      "NJA00412", "NJA00404", "NJA00405", "NJA00406", "NJA00408", "NJA00409",
      "NJA00413", "NJA00602", "NJA00608", "NJA00609", "NJA00801", "NJA01202",
      "NJA01203", "NJA01205", "NJA01301", "NJA01302", "NJA01303", "NJA01401"
    ),
    NA, GFP)
  ) %>%
  annotateMouse(
    .,
    plate = c(
      "NJA00110", "NJA00101", "NJA00111", "NJA00107", "NJA00102", "NJA00103",
      "NJA00104", "NJA00109", "NJA00201", "NJA00204", "NJA00205", "NJA00206",
      "NJA00402", "NJA00403", "NJA00411", "NJA00412", "NJA00404", "NJA00405",
      "NJA00406", "NJA00408", "NJA00409", "NJA00413", "NJA00602", "NJA00608",
      "NJA00609", "NJA00801", "NJA01202", "NJA01203", "NJA01205", "NJA01301",
      "NJA01302", "NJA01303", "NJA01401"
    ),
    mouse = c(
      1, 1, 1, 1, 1, 1,
      1, 1, 2, 2, 2, 2,
      4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 6, 6,
      6, 8, 12, 12, 12,
      13, 13, 13, 14
    )
  ) %>%
  annotateTissue(
    .,
    plate = c(
      "NJA00101", "NJA00102", "NJA00103", "NJA00104", "NJA00107", "NJA00109",
      "NJA00110", "NJA00111", "NJA00201", "NJA00204", "NJA00205", "NJA00206",
      "NJA00402", "NJA00403", "NJA00411", "NJA00412", "NJA00404", "NJA00405",
      "NJA00406", "NJA00408", "NJA00409", "NJA00413", "NJA00602", "NJA00608",
      "NJA00609", "NJA00801", "NJA01202", "NJA01203", "NJA01205", "NJA01301",
      "NJA01302", "NJA01303", "NJA01401"
    ),
    tissue = c(
      "SI", "SI", "SI", "SI", "SI", "colon",
      "colon", "colon", "colon", "colon", "colon", "colon",
      "SI", "SI", "colon", "colon", "SI", "SI", "SI",
      "colon", "colon", "colon", "SI", "colon",
      "colon", "SI", "SI", "colon", "colon", "SI",
      "SI", "colon", "colon"
    )
  ) %>%
  dplyr::mutate(cellNumber = dplyr::if_else(
    stringr::str_detect(sample, "^s"), "Singlet", "Multiplet")
  ) %>%
  dplyr::mutate(filtered = dplyr::if_else(sample %in% colnames(counts), FALSE, TRUE))

#rename and save
countsMgfp <- counts
countsMgfpERCC <- countsERCC
countsMgfpMeta <- plateData

if(all(!colnames(countsMgfp) %in% countsMgfpMeta$sample)) {
  stop("all counts data not present in meta data")
}
if(all(!colnames(countsMgfpERCC) %in% countsMgfpMeta$sample)) {
  stop("all ercc data not present in meta data")
}

save(
  countsMgfp,
  countsMgfpERCC,
  countsMgfpMeta,
  file = "./data/countsMgfp.rda",
  compress = "bzip2"
)
