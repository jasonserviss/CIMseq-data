#run from package root
#source('./inst/rawData/countsSorted2/Sorted.multiplets.R')

packages <- c("sp.scRNAseqData", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

projectName <- "Sorted.cell.line.multiplets_SCM"
shortName <- "SCM"
cat(paste0('Processing ', projectName, '\n'))

googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")
Meta <- getMetadata(projectName)
if("Missing" %in% colnames(Meta)) {
  Meta <- Meta %>%
    filter(is.na(Missing) | Missing == FALSE) %>%
    select(-Missing)
}

countData <- getCountsData(projectName)

#move genes to rownames
countData <- moveGenesToRownames(countData)

#remove .htseq suffix
if(any(grepl("htseq", colnames(countData)))) countData <- removeHTSEQsuffix(countData)

#label singlets and multiplets ids should include SINGLET samples only
singlets <- Meta[Meta$cellNumber == "Singlet", "sample"][[1]]
singlets <- gsub("^..(.*)", "\\1", singlets)
counts <- labelSingletsAndMultiplets(countData, singlets)

#extract ERCC
ercc <- detectERCCreads(counts)
countsERCC <- counts[ercc, ]
counts <- counts[!ercc, ]

#remove non-genes
counts <- counts[!detectNonGenes(counts), ]
namesPreFilter <- colnames(counts)

data <- filterCountsData(
  counts, countsERCC, geneMinCount = 0, cellMinCount = 4e4, geneName = "ACTB",
  quantileCut.hk = 0.01, quantileCut.ercc = 0.95
)

#add filtered column to Meta
Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
  sample %in% colnames(data[[1]]),
  FALSE, TRUE
))

#check all count samples in meta and vice versa
c1 <- all(!Meta$sample %in% namesPreFilter)
c2 <- all(!colnames(data[[1]]) %in% Meta$sample)
if(c1 & c2) {
  stop("all counts data not present in meta data")
}

#rename
Counts <- data[[1]]
CountsERCC <- data[[2]]

#save .rda
saveRDA(projectName, Counts, CountsERCC, Meta)

#upload .txt
processedDataUpload(projectName, Counts, CountsERCC, Meta)
