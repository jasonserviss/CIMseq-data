#run from package root
#source('inst/rawData/countsMgfp/count_table_180112.R')

#load counts
#remove columns with "SI1.Singlet.1" in colnames. Early miSeq. Duplicated samples.
#several samples are only NA. find with apply(counts, 2, function(x) all(is.na(x)))
#Regev data identified by "SRR" in colnames

path <- './inst/rawData/countsMgfp/count_table_180112.txt'
counts <- read.table(path, header = TRUE, sep = "\t")
bool1 <- !grepl("SI1.Singlet.1", colnames(counts)) #old miSeq data
bool2 <- !grepl("SRR", colnames(counts)) #Regev data
bool3 <- !apply(counts, 2, function(x) all(is.na(x))) #samples without counts
bool4 <- "HGN" %in% colnames(counts)
counts <- counts[, bool1 & bool2 & bool3 & bool4]

#check for NAs
if(sum(is.na(counts)) > 0) {
  stop("NAs in counts data.")
}

#move genes to rownames
counts <- moveGenesToRownames(counts)

#remove "htseq" suffix
counts <- removeHTSEQsuffix(counts)

#label singlets and multiplets
ids <- c(
  "Singlet", "NJA00102", "NJA00103", "NJA00104",
  "NJA00109", "NJA00204", "NJA00205", "NJA00206"
)
counts <- labelSingletsAndMultiplets(counts, ids)

#remove old singlet multiplet specification
colnames(counts) <- gsub("Singlet\\.", "", colnames(counts))
colnames(counts) <- gsub("Doublet\\.", "", colnames(counts))

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts), ]

#remove low quality cells
lqc <- detectLowQualityCells(counts, geneName = "Actb")
counts <- counts[, lqc]
countsERCC <- countsERCC[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)
countsERCC <- convertCountsToMatrix(countsERCC)

#rename and save
countsMgfp <- counts
countsMgfpERCC <- countsERCC
save(
  countsMgfp,
  countsMgfpERCC,
  file = "./data/countsMgfp.rda",
  compress = "bzip2"
)
