#run from package root
#source('inst/rawData/fetalPancreas/fetalPancreasCounts.R')

library(sp.scRNAseqData)
cat('Processing fetalPancreas.\n')

#fetalPancreasCounts
#Should be 131 singlets and 69 multiplets.
#Note: I never got the raw unfiltered counts.txt file for this from Martin.

googledrive::drive_auth(oauth_token = "data/gd.rds")

#download data
googledrive::drive_download(
  file = 'fetalPancreasCounts.txt',
  path = './inst/rawData/fetalPancreas/fetalPancreasCounts.txt',
  overwrite = TRUE
)

#read counts data
counts <- read.table(
  'inst/rawData/fetalPancreas/fetalPancreasCounts.txt',
  header = TRUE,
  sep = "\t"
)

#add prefix indicating singlets and multiplets
counts <- labelSingletsAndMultiplets(counts, "1000102901")

#extract ercc reads
ercc <- detectERCCreads(counts)
countsErccFp <- counts[ercc, ]
countsFp <- counts[!ercc, ]

#save data
save(
  countsFp,
  countsErccFp,
  file = 'data/fetalPancreasCounts.rda',
  compress = 'bzip2'
)

