
#run method
packages <- c(
  "sp.scRNAseq",
  "sp.scRNAseqTesting",
  "tidyverse"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

sng <- str_detect(colnames(countsSorted2), "^s")

#create counts objects
cObjSng <- spCounts(countsSorted2[, sng], countsSortedERCC2[, sng])
cObjMul <- spCounts(countsSorted2[, !sng], countsSortedERCC2[, !sng])

#spUnsupervised
uObj <- spUnsupervised(cObjSng)

#rename classes
positions <- str_extract(colnames(getData(cObjSng, "counts")), "...$")
newClass <- case_when(
  positions == "E03" ~ "A375", #what is this? sorting issue?
  positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("01", "02", "03", "04")) ~ "HOS",
  positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("05", "06", "07", "08")) ~ "HCT116",
  positions %in% paste0(sort(rep(LETTERS[1:8], 4)), c("09", "10", "11", "12")) ~ "A375",
  TRUE ~ "error"
)
corresp <- getData(uObj, "classification") %>%
tibble(oldClass = ., newClass = newClass) %>%
distinct()

classification(uObj) <- newClass

gm <- getData(uObj, "groupMeans")
colnames(gm) <- pull(corresp, newClass)[match(colnames(gm), pull(corresp, oldClass))]
groupMeans(uObj) <- gm

tm <- getData(uObj, "tsneMeans")
tm$classification <- pull(corresp, newClass)[match(tm$classification, pull(corresp, oldClass))]
tsneMeans(uObj) <- tm

#permute and calculate class means
groupMeanPermutations <- permuteMeans(cObjMul, uObj)

save(
  groupMeanPermutations,
  file = './data/groupMeanPermutations.rda',
  compress = "bzip2"
)

#source('./inst/rawData/distPermutations.R')
