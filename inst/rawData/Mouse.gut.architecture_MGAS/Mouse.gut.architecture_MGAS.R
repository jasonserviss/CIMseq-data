#run from package root
#source('inst/rawData/countsMgfp/Mouse.gut.architecture.R')

packages <- c("sp.scRNAseqData", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

projectName <- "Mouse.gut.architecture.segments_MGAS"
shortName <- "MGAS"
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

#label singlets and multiplets ids should include SINGLET samples only
colnames(countData) <- paste("m", colnames(countData), sep = ".")

#extract ERCC
ercc <- detectERCCreads(countData)
CountsERCC <- countData[ercc, ]
Counts <- countData[!ercc, ]

#remove non-genes
Counts <- Counts[!detectNonGenes(Counts), ]

#filter counts
data <- filterCountsData(
  Counts, CountsERCC, geneMinCount = 0, cellMinCount = 1.8e4, geneName = "Actb",
  quantileCut.hk = 0.01, quantileCut.ercc = 0.99
)

#add filtered column to Meta
Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
  sample %in% colnames(data[[1]]),
  FALSE, TRUE
))

#check all count samples in meta and vice versa
c1 <- all(!Meta$sample %in% colnames(Counts))
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
