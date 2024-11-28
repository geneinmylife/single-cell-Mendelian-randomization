library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(Matrix)
library(stringr)
library(ggsignif)
library(ggpubr)
library(viridis)
library(data.table)


file_list <- c('H1','H2','N1','N2','N3','T1','T2')
data_list <- list()
for(i in 1:length(file_list)){
  tmp_data <- myRead10X(data.dir = 'path',prefixed = file_list[i])
  colnames(tmp_data) <- paste(file_list[i],colnames(tmp_data),sep = '_')
  data_list[[i]] <- tmp_data
}
colnames(data_list[[1]])

data_total <- do.call(cbind,data_list)
obj1 <- CreateSeuratObject(counts = data_total,min.cells = 100)
obj1$percent.mt <- PercentageFeatureSet(obj1,pattern = "^MT-")
obj1$disease <- str_split_fixed(colnames(obj1),'_',2)[,1]
obj1 <- subset(obj1,percent.mt<20)
obj1 <- SCTransform(obj1)
obj1 <- RunPCA(obj1)
obj1 <- RunHarmony(obj1,group.by.vars='disease')
obj1 <- RunUMAP(obj1,dims = 1:15,reduction = 'harmony')
DimPlot(obj1)

obj1 <- FindNeighbors(obj1,dims = 1:15,reduction = 'harmony')
obj1 <- FindClusters(obj1,resolution = 0.1)
DimPlot(obj1)
FeaturePlot(obj1,features = c('CD3D','CD4','CD8A'))


obj3 <- subset(obj1,seurat_clusters==0)
obj3 <- RunPCA(obj3)
obj3 <- RunHarmony(obj3,group.by.vars='disease')
obj3 <- RunUMAP(obj3,dims = 1:15,reduction = 'harmony')
obj3 <- FindNeighbors(obj3,dims = 1:15,reduction = 'harmony')
obj3 <- FindClusters(obj3,resolution = 0.1)
DimPlot(obj3)
FeaturePlot(obj3,features = c('CD4','CD8A','FOXP3'))
obj3$disease2 <- 'Tumor'
obj3$disease2[grepl('N',obj3$disease)] <- 'Normal'
obj3@active.ident <- factor(obj3$disease2)


Tconv_obj <- subset(obj3,seurat_clusters==0)
Treg_obj <- subset(obj3,seurat_clusters==3)

Treg_marker1 <- FindMarkers(Treg_obj,ident.1 = 'Tumor',min.pct = 0,logfc.threshold = 0)
conv_marker1 <- FindMarkers(Tconv_obj,ident.1 = 'Tumor',min.pct = 0,logfc.threshold = 0)

