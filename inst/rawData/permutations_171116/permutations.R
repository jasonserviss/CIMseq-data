#run from package root
#source('inst/rawData/permutations_171116/permutations.R')

packages <- c('sp.scRNAseq', 'sp.scRNAseqTesting', 'sp.scRNAseqData', 'stringr', 'dplyr', 'tibble')
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#run the method
sng <- str_detect(colnames(countsSorted2), "^s")
cObjSng <- spCounts(countsSorted2[, sng], countsSortedERCC2[, sng])
cObjMul <- spCounts(countsSorted2[, !sng], countsSortedERCC2[, !sng])
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

#run real
cellNrs <- estimateCells(cObjSng, cObjMul)
sObjPermutations <- spSwarm(spCounts = cObjMul, spUnsupervised = uObj, distFun = "dtsnCellNum", maxiter = 10, swarmsize = 50, cores = 12, norm = TRUE, cellNumbers = cellNrs, e = 0.0025)

#run permutations
permutations <- permuteSwarm(spCountsMul = cObjMul, spUnsupervised = uObj, distFun = sp.scRNAseq:::dtsnCellNum, maxiter = 10, swarmsize = 50, cores = 12, norm = TRUE, iter = 10000, cellNumbers = cellNrs, e = 0.0025)

save(sObjPermutations, permutations, file = "data/permutations.rda", compress = "bzip2")


