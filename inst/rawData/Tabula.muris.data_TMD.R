#run from package root
#source('./inst/rawData/countsRegev/Regev.small.intestine.R')
#GSE92332

packages <- c("CIMseq.data", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

TMD <- function(upload = TRUE, save = TRUE) {
  projectName <- "Tabula.muris.data_TMD"
  shortName <- "TMD"
  cat(paste0('Processing ', projectName, '\n'))
  
  googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")
  Meta <- getMetadata(projectName)
  if("Missing" %in% colnames(Meta)) {
    Meta <- Meta %>%
      filter(is.na(Missing) | Missing == FALSE) %>%
      select(-Missing)
  }
  
  countData <- getCountsData(projectName)
  
  #move gene names to rownames
  countData <- moveGenesToRownames(countData)
  
  #remove .htseq suffix
  #countData <- removeHTSEQsuffix(countData)
  
  #annotate all samples as singlets
  colnames(countData) <- paste0("s.", colnames(countData))
  
  #extract ERCC reads
  ercc <- detectERCCreads(countData)
  countsERCC <- countData[ercc, ]
  counts <- countData[!ercc, ]
  
  #remove non-genes
  counts <- countData[!detectNonGenes(counts), ]
  
  #filter bad genes and cells
  data <- filterCountsData(
    counts, countsERCC, geneMinCount = 0, cellMinCount = 1e4, geneName = "Actb",
    quantileCut.hk = 0.01, quantileCut.ercc = 0.99
  )
  
  #add filtered column to Meta
  Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
    sample %in% colnames(data[[1]]),
    FALSE, TRUE
  ))
  
  #rename
  Counts <- data[[1]]
  CountsERCC <- data[[2]]
  
  #save .rda
  if(save) saveRDA(projectName, Counts, CountsERCC, Meta)
  
  #upload .txt
  if(upload) processedDataUpload(projectName, Counts, CountsERCC, Meta)
}
