#run from package root
#script('inst/rawData/fetalPancreas/fetalPancreasCounts.R')

library(sp.scRNAseqData)

#fetalPancreasCounts
#Should be 131 singlets and 69 multiplets.
#Note: I never got the raw unfiltered counts.txt file for this from Martin.

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

