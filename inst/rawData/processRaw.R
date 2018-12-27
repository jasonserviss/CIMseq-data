#source('./inst/rawData/processRaw.R')

#run scripts to generate data/.rda files
ignore <- c(
  './inst/rawData/permutations_171116/permutations.R',
  './inst/rawData/syntheticData/syntheticData.R',
  './inst/rawData/countsSorted1/countsSorted1_171018.R',
  './inst/rawData/Human.fetal.pancreas_HFP.R',
  './inst/rawData/Aom.dss.mouse_ADM.R',
  './inst/rawData/processRaw.R'
)

processRaw <- function(
  basePath = './inst/rawData', ignore = NULL, upload = TRUE, save = TRUE
){
  paths <- list.files(
    path = basePath, recursive = TRUE, pattern = '\\.R$', full.names = TRUE
  )

  if(!is.null(ignore)) paths <- paths[!paths %in% ignore]
  funNames <- gsub(".*_(.*)\\.R", "\\1", basename(paths))

  trash <- purrr::map2(paths, funNames, function(p, fn) {
    source(p)
    c.fn <- get(fn)
    c.fn(upload = upload, save = save)
  })
}
