
#A function to get geneset scores of scRNA-Seq data with AUCell Package
#https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#exploring-the-cell-assignment-table-heatmap 


BiocManager::install("AUCell")
library(AUCell)
library(Seurat)


#The function takes 2 argument
#1.Seurat object
#2.Geneset vector



func_au<- function (seu_obj, gene_list) {
  #get counts data from the seurat object
  seu_counts<- GetAssayData(object = seu_obj, slot = "counts")
  #build cell rankings
  cells_rankings <- AUCell_buildRankings(seu_counts,splitByBlocks=TRUE)
  #calculate AUCell score
  cells_AUC <- AUCell_calcAUC(gene_list,
                              cells_rankings,aucMaxRank=nrow(cells_rankings)*0.05)
  #plot the AUC histogram
  cell_assignments <- AUCell_exploreThresholds(cells_AUC,plotHist=TRUE,assign=TRUE)
  #get the AUCell score 
  #tranpose, since we want the cells on the rows
  AUCscore_mtx<-t(getAUC(cells_AUC))
  return(AUCscore_mtx)
  
 }
