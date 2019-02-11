#run from package root
#source('inst/rawData/countsMgfp/Mouse.gut.architecture.R')

packages <- c("CIMseq.data", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

MGA <- function(upload = TRUE, save = TRUE) {
  projectName <- "Mouse.gut.architecture_MGA"
  shortName <- "MGA"
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
  sng <- grepl("^s", colnames(Counts))
  data.s <- filterCountsData(
    Counts[, sng], CountsERCC[, sng], geneMinCount = 0, cellMinCount = 1.8e4, 
    geneName = "Actb", quantileCut.hk = 0.01, quantileCut.ercc = 0.99
  )
  data.m <- filterCountsData(
    Counts[, !sng], CountsERCC[, !sng], geneMinCount = 0, cellMinCount = 2.3e4, 
    geneName = "Actb", quantileCut.hk = 0.01, quantileCut.ercc = 0.99
  )
  genes <- intersect(rownames(data.s[[1]]), rownames(data.m[[1]]))
  samples <- c(colnames(data.s[[1]]), colnames(data.m[[1]]))
  
  #add filtered column to Meta
  Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
    sample %in% samples, FALSE, TRUE
  ))

  #check all count samples in meta and vice versa
  c1 <- all(!Meta$sample %in% colnames(Counts))
  c2 <- all(!samples %in% Meta$sample)
  if(c1 & c2) {
    stop("all counts data not present in meta data")
  }

  #rename
  Counts <- cbind(data.s[[1]][genes, ], data.m[[1]][genes, ])
  CountsERCC <- cbind(data.s[[2]], data.m[[2]])

  #save .rda
  if(save) saveRDA(projectName, Counts, CountsERCC, Meta)

  #upload .txt
  if(upload) processedDataUpload(projectName, Counts, CountsERCC, Meta)
}
