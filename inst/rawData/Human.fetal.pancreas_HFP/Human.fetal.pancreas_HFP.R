#run from package root
#source('inst/rawData/fetalPancreas/fetalPancreasCounts.R')

packages <- c("sp.scRNAseqData", "EngeMetadata", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

HFP <- function(upload = TRUE, save = TRUE) {
  projectName <- "Human.fetal.pancreas_HFP"
  shortName <- "HFP"
  cat(paste0('Processing ', projectName, '\n'))

  #fetalPancreasCounts
  #Should be 131 singlets and 69 multiplets.
  #Note: I never got the raw unfiltered counts.txt file for this from Martin.

  googledrive::drive_auth(oauth_token = "inst/extData/gd.rds")

  Meta <- getMetadata(projectName)
  if("Missing" %in% colnames(Meta)) {
    Meta <- Meta %>%
      filter(is.na(Missing) | Missing == FALSE) %>%
      select(-Missing)
  }

  #read counts data
  countData <- getCountsData(projectName)

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

  #add filtered column to Meta
  Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
    sample %in% colnames(Counts),
    FALSE, TRUE
  ))

  #check all count samples in meta and vice versa
  c1 <- all(!Meta$sample %in% colnames(Counts))
  c2 <- all(!colnames(Counts) %in% Meta$sample)
  if(c1 & c2) {
    stop("all counts data not present in meta data")
  }

  #save .rda
  if(save) saveRDA(projectName, Counts, CountsERCC, Meta)

  #upload .txt
  if(upload) processedDataUpload(projectName, Counts, CountsERCC, Meta)
}
