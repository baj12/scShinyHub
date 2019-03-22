
# smaller data set with some random data
load("Examples/PBMC-Apheresis.new.Rds")

library(SingleCellExperiment)
require(irlba)

class(gbm)
gbmMat = as.matrix(gbm)
object.size(gbm)
object.size(gbmMat)
counts <- as(gbmMat, "dgCMatrix")

scEx = SingleCellExperiment(assay = list(counts=counts),
                            colData = pData(gbm),
                            rowData = fData(gbm))
class(scEx)
object.size(scEx)
class(assays(scEx)[[1]])

system.time({scEx = scater::runPCA(scEx,ncomponents = 10, method = "irlba",
               ntop = 500, exprs_values = "counts")})
pca = reducedDim(scEx, "PCA")
attr(pca,"percentVar")

system.time({irlaPCA = irlba::prcomp_irlba(assays(scEx)[[1]], n = 10, 
                         retx = TRUE, 
                         center = Matrix::colMeans(assays(scEx)[[1]]),
                         fastpath=FALSE,
                         scale. = FALSE)})
cent = Matrix::colMeans(assays(scEx)[[1]])
rv <- rowVars(assays(scEx)[[1]])
system.time({prPCA = prcomp_irlba(assays(scEx)[[1]], n = 10, 
                         retx = TRUE, 
                         center = cent,
                         fastpath=TRUE,
                         scale. = FALSE)})

save(file = "scEx.Rds", list = c("scEx"))

summary(scEx)
colData(scEx)
rownames(scEx)
names(scEx)
colnames(scEx)
exprs(scEx) # doesn't work
sum(is.infinite(scEx))
assaysSCEX = (assays(scEx))
assays(scEx)[[1]]
sampleIdx = sample.int(ncol(gbm), 1000)
Matrix::colSums(scEx > 0, na.rm = TRUE)
ncol(countES) <= 1 | nrow(countES) < 1
colnames(dataTables$scEx)
selCols <- scEx[ids, ] > 0
pData(countES)$sampleNames
gbm_bcnorm <- (newGeneBCMatrix(A, fData(scEx), pData(scEx),
                               template = scEx
))
colData(scEx)$sampleNames
colData(gbmNew) <- pD








gbm = gbm[,sampleIdx]

pd = pData(gbm)
rownames(pd)[1:500] = gsub("-1", "-samp1", rownames(pd)[1:500])
rownames(pd)[501:1000] = gsub("-1", "-samp2", rownames(pd)[501:1000])
pd$sampleNames = as.character(pd$sampleNames)
pd$sampleNames[1:500] = "samp1"
pd$sampleNames[501:1000] = "samp2"
pd$sampleNames = factor(pd$sampleNames)
pd$barcode = rownames(pd)

pd$randomNr = runif(nrow(pd),0,10)
pd$randomNr2 = runif(nrow(pd),0,10000)

pData(gbm) = pd
ex = exprs(gbm)
colnames(ex) = rownames(pd)
# exprs(gbm) = ex
gbm <- newGeneBCMatrix(mat = as(ex, "dgTMatrix"), fd = fData(gbm), pd = pd)

