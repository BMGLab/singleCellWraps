# How can ı merge my seurat objects?
library(Seurat)
library(tidyverse)
#first object
seub1 <- CreateSeuratObject(counts = cnt_b1, project = "B01", min.cells = 3, min.features = 200)
seub1@meta.data$Outcome <- c("Responder") #you can easily add any feature on metadata with that way
rm(cnt_b1) #optional, you should remove your existing objects for avoid high cpu usage
gc()
#second object
seub2 <- CreateSeuratObject(counts = cnt_b2, project = "B02", min.cells = 3, min.features = 200)
seub2@meta.data$Outcome <- c("Non-Responder")
head(seub2@meta.data)
rm(cnt_b2)
gc()

#merging
seu <- merge(seub1, y = seub3, add.cell.ids = c("B01", "B03"), project = "DLIstudy")

#standard seurat pipeline
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:10)
DimPlot(seu, reduction = "umap", group.by = "Outcome")
DimPlot(seu, reduction = "pca", group.by = "Outcome")
