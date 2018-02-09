#run from package root
#source('./inst/rawData/countsRegev/count_table_180112.R')
#GSE92332
#load counts
#Subset Regev data using SRR in colnames
path <- './inst/rawData/countsMgfp/count_table_180112.txt'
counts <- read.table(path, header = TRUE, sep = "\t")
bool1 <- grepl("SRR", colnames(counts))
bool2 <- colnames(counts) == "HGN"
counts <- counts[, bool1 | bool2]

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
#there are no ercc in the data
table(detectERCCreads(counts))

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]

#remove low quality genes
counts <- counts[detectLowQualityGenes(counts), ]

#remove low quality cells
lqc <- detectLowQualityCells(
  counts,
  geneName = "Actb",
  mincount = 5e4,
  quantileCut = 0.01
)
counts <- counts[, lqc]

#coerce to matrix
counts <- convertCountsToMatrix(counts)

#rename and save
countsRegev <- counts

save(
  countsRegev,
  file = "./data/countsRegev.rda",
  compress = "bzip2"
)
