#run from package root
#source('inst/rawData/countsMgfp/Mouse.gut.architecture.R')

packages <- c("sp.scRNAseqData", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

projectName <- "Mouse.gut.architecture"
cat(paste0('Processing ', projectName, '\n'))

googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")
Meta <- getMetadata(projectName)
if("Missing" %in% colnames(Meta)) {
  Meta <- Meta %>%
    filter(is.na(Missing) | Missing == FALSE) %>%
    select(-Missing)
}

countData <- getCountsData(projectName)
colnames(countData) <- renameMgfpSamples(colnames(countData))

#move genes to rownames
countData <- moveGenesToRownames(countData)

#label singlets and multiplets ids should include SINGLET samples only
singlets <- Meta[Meta$cellNumber == "Singlet", "sample"][[1]]
singlets <- gsub("^..(.*)", "\\1", singlets)
countData <- labelSingletsAndMultiplets(countData, singlets)

#extract ERCC
ercc <- detectERCCreads(countData)
CountsERCC <- countData[ercc, ]
Counts <- countData[!ercc, ]

#remove non-genes
Counts <- Counts[!detectNonGenes(Counts), ]

#filter counts
data <- filterCountsData(
  Counts, CountsERCC, geneMinCount = 0, cellMinCount = 1e4, geneName = "Actb",
  quantileCut = 0.01, percentile = 0.99
)

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
