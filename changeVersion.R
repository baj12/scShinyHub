require(Biobase)
library(SingleCellExperiment)
require(irlba)

infile = "~/Downloads/AVC-mm.Rds"

inFiles = dir(path = ".", pattern = "*.Rds$", full.names = T, recursive = T)
inFiles = dir(path = "~/GoogleDrive/", pattern = "*.Rds$", full.names = T, recursive = T)

inFiles = dir(path = "~/Rstudio/", pattern = "*.Rds$", full.names = T, recursive = T)
infile = inFiles[1]
for (infile in inFiles) {
  if (tools::file_ext(tools::file_path_sans_ext(infile)) == "v2") next()
  outfile = paste0(tools::file_path_sans_ext(infile) , ".v2.Rds")
  # if (file.exists(outfile)) next()
  cat(file = stderr(), paste(infile, "\n"))
  
  tryCatch({rm(list = c("gbm", "featuredata"))},
           error = function(e){},
           warning = function(e){})
  base::load(infile)
  if (!"gbm" %in% ls()) next()
  pd = pData(gbm)
  fd = fData(gbm)
  # cat(file = stderr(), paste(featuredata[which(!rownames(featuredata) %in% rownames(fd)),], collapse = "\n"))
  fd = cbind(fd, featuredata[rownames(fd),which(!colnames(featuredata) %in% c("Chromosome.Name", "Gene.Start..bp.", "Gene.End..bp." ) )])
  pd$sampleNames = as.factor(pd$sampleNames)
  counts <- as(exprs(gbm), "dgCMatrix")
  scEx = SingleCellExperiment(assay = list(counts=counts),
                              colData = pd,
                              rowData = fd)
  
  save(file = outfile, list = c("scEx"))
}
