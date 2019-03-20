library(Seurat)
library(dplyr)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size

sparse.size <- object.size(x = pbmc.data)
sparse.size


dense.size/sparse.size


pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")
object.size(pbmc)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)


pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)


pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))



pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)



PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)




VizPCA(object = pbmc, pcs.use = 1:2)


PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)



pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)



PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)



PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)



pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)



JackStrawPlot(object = pbmc, PCs = 1:12)



PCElbowPlot(object = pbmc)


pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)


PrintFindClustersParams(object = pbmc)


pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)



TSNEPlot(object = pbmc)

object.size(pbmc)
saveRDS(pbmc, file = "pbmc_tutorial.rds")
pbmc@raw.data
pbmc@meta.data # some of the pData
names(pbmc@dr) # projections / dimension reductions
pbmc@hvg.info # highly variable genes ( could be a group of genes?)
pbmc@imputed
pbmc@cluster.tree
pbmc@snn
pbmc@kmeans
pbmc@spatial
pbmc@misc
