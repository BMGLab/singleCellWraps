#Essential Seurat Workflow

#Seurat workflow takes raw single-cell expression data and aims to find clusters within the data

#You must have count matrix, features and cell barcodes present in the folder:
pbmc.counts <- Read10X(data.dir = "~/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Create the Seurat object with the raw data:
pbmc <- CreateSeuratObject(counts = pbmc.counts)
#QC to filter out doublets and the dead cells:
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#Normalization: to normalize the feature expression measurements for each cell by the total expression, 
#multiplies this by a scale factor (10,000 by default), and log-transforms the result.
pbmc <- NormalizeData(object = pbmc)
#Variable Features:highly expressed in some cells, and lowly expressed in others
#These genes helps to highlight biological signal in single-cell dataset in downstream analysis 
pbmc <- FindVariableFeatures(object = pbmc)
#ScaleData: linear transformation which is a standard pre-processing step prior to dimensional reduction like PCA
pbmc <- ScaleData(object = pbmc)
#PCA
pbmc <- RunPCA(object = pbmc)
#To select which PCs will be used for the analysis
ElbowPlot(pbmc)
#clustering: K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, 
#and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities.
#Then apply modularity optimization techniques like Louvain algorithm 
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
#To visualize and explore these dataset
#Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots
pbmc <- RunUMAP(object = pbmc)
DimPlot(object = pbmc, reduction = "umap")