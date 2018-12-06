# convert an old zip file to an Rds file that can be read
rm(list=ls())
zipDir = "~/Desktop/tartour/report.sected/"
rdsPath = "/Users/bernd/GoogleDrive/pasteur/Tartour/"
outFile = paste0(zipDir, "sessionData.Rds")

library(readr)

normalizedCounts <- read_csv(paste0(zipDir,"normalizedCounts.csv"))
r <- as.data.frame(normalizedCounts)
rownames(r) <- r[, 1]
normalizedCounts <- r[, -1]
gqcProjections <- read_csv(paste0(zipDir,"gqcProjections.csv"))
r <- as.data.frame(gqcProjections)
rownames(r) <- r[, 1]
gqcProjections <- r[, -1]
load(paste0(zipDir, "sessionData.RData"))
# we need to know where to find the input file:
load(paste0(rdsPath, params$file1$name))

cells = rownames(gqcProjections)

genes = rownames(normalizedCounts)
mat <- exprs(gbm)[genes, cells]
fd <- fData(gbm)[genes,]
pd <- pData(gbm)[cells,]
pd$sampleNames = factor(pd$sampleNames)
gbm = newGeneBCMatrix(mat, fd, pd)
colnames(fData(gbm))
table(pd$sampleNames)

sum(c("id", "symbol") %in% colnames(fData(gbm)))
save(file = outFile, list = c("gbm", "featuredata"))
