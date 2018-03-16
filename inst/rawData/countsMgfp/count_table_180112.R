#run from package root
#source('inst/rawData/countsMgfp/count_table_180112.R')

library(sp.scRNAseqData)
library(stringr)
library(dplyr)

#load counts
#several samples are only NA. find with table(apply(counts, 2, function(x) all(is.na(x))))
#Regev data identified by "SRR" in colnames

#path <- './inst/rawData/countsMgfp/count_table_180112.txt'
path <- './inst/rawData/countsMgfp/counts_table_180316.txt'
counts <- read.table(path, header = TRUE, sep = "\t")
#bool2 <- !grepl("SRR", colnames(counts)) #Regev data
#bool4 <- "HGN" %in% colnames(counts)
#counts <- counts[, bool2 & bool4]
#bool1 <- !grepl("SI1.Singlet.1", colnames(counts)) #old miSeq data (should be removed now)
#bool3 <- !apply(counts, 2, function(x) all(is.na(x))) #samples without counts (should be resolved now)
#counts <- counts[, bool1 & bool3]

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#remove "htseq" suffix
#counts <- removeHTSEQsuffix(counts)

#label singlets and multiplets
ids <- c(
  "Singlet", "NJA00102", "NJA00103", "NJA00104",
  "NJA00109", "NJA00204", "NJA00205", "NJA00206"
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
counts <- counts[detectLowQualityGenes(counts), ]

#remove low quality cells
lqc <- detectLowQualityCells(
  counts,
  geneName = "Actb",
  mincount = 1e5,
  quantileCut = 0.01
)
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#prepare metadata
plateData <- loadMetaData('./inst/rawData/countsMgfp/counts_Mgfp_meta.txt') %>%
mutate(sample = removeHTSEQsuffix(sample)) %>%
mutate(sample = labelSingletsAndMultiplets(
  sample,
  c(
    "Singlet", "NJA00102", "NJA00103", "NJA00104",
    "NJA00109", "NJA00204", "NJA00205", "NJA00206"
))) %>%
mutate(sample = renameMgfpSamples(sample)) %>%
annotatePlate(.) %>%
annotateRow(.) %>%
annotateColumn(.) %>%
annotateGFP(
  .,
  plate = list(c("NJA00111"), c("NJA00110"), c("NJA00103"), c("NJA00104"), c("NJA00201"), c("NJA00110")),
row = list(1:8, 1:8, 1:8, 1:8, 1:8, 1:8),
column = list(1:12, 1:12, 1:6, 1:6, 1:12, 7:12)
) %>%
mutate(GFP = if_else(plate %in% c("NJA00204", "NJA00205", "NJA00206"), NA, GFP)) %>%
annotateMouse(
  .,
  plate = c(
    "NJA00110", "NJA00101", "NJA00111", "NJA00107", "NJA00102", "NJA00103",
    "NJA00104", "NJA00109", "NJA00201", "NJA00204", "NJA00205", "NJA00206"
  ),
  mouse = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2)
) %>%
annotateTissue(
  .,
  plate = c(
    "NJA00110", "NJA00101", "NJA00111", "NJA00107", "NJA00102", "NJA00103",
    "NJA00104", "NJA00109", "NJA00201", "NJA00204", "NJA00205", "NJA00206"
  ),
  tissue = c(
    "colon", "SI", "colon", "colon", "SI", "SI", "SI", "SI", "colon", "colon",
    "colon", "colon", "colon"
  )
) %>%
mutate(cellNumber = if_else(str_detect(sample, "^s"), "Singlet", "Multiplet"))

#rename and save
countsMgfp <- counts
countsMgfpERCC <- countsERCC
countsMgfpMeta <- plateData
save(
  countsMgfp,
  countsMgfpERCC,
  countsMgfpMeta,
  file = "./data/countsMgfp.rda",
  compress = "bzip2"
)
