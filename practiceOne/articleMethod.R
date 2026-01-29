#This script uses the cell clustering method described in this article: https://www.nature.com/articles/s41593-021-00938-x#Sec43
install.packages("Seurat")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("harmony")
install.packages('leidenbase')
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
#Required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
#raw data
raw_data
#seed
set.seed(67)

data<- raw_data
data <- NormalizeData(data, normalization.method = "LogNormalize")
data <- ScaleData(data, features = rownames(data))
seurat_obj <- RunPCA(data, features = rownames(data), verbose = FALSE)

#Elbow plot 1
ElbowPlot(seurat_obj, ndims = 50)

#Figure 1
#Elbow plot 2, helps with more clarity
ElbowPlot(seurat_obj, ndims = 20)

#Decided to use 10 dimensions

rdOneDims <- 10

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:rdOneDims, k.param = 25)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:rdOneDims)
#UMAP generated

#Figure 2
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

#Markers for different cell types that were laid out in the paper
round1Markers <- c("Ppp1r1b",  # MSN
                  "Resp18",    # Interneuron
                  "Snap25",    # Other neuronal
                  "Aldoc",     # Astrocytes
                  "C1qa",      # Microglia
                  "Mobp",      # Oligodendrocytes
                  "Pdgfra",    # OPCs
                  "Rgs5",      # Endothelial
                  "Ccdc153")   # Ependymal

FeaturePlot(seurat_obj, features = round1Markers, ncol = 3)

#Figure 3
DotPlot(seurat_obj, features = round1Markers) + RotatedAxis()

#Figure 4
#Microglia
DotPlot(seurat_obj, features = "C1qa")
seurat_obj <- RenameIdents(seurat_obj, '10' = "Microglia")

#Figure 5
#Endothelial
DotPlot(seurat_obj, features = "Rgs5")
seurat_obj <- RenameIdents(seurat_obj, '14' = "Endothelial")

#Figure 6
#Oligodendrocyte Precursors
DotPlot(seurat_obj, features = "Pdgfra")
seurat_obj <- RenameIdents(seurat_obj, '13' = "Oligodendrocyte Precursors")

feats <- c("Drd1", "Drd2")
DotPlot(seurat_obj, features = feats)
FeaturePlot(seurat_obj, features = feats)


FeaturePlot(seurat_obj, features = "C1qa")

save.image(file = "articleMethod.RData")

