---
title: "ScRNAseq_training"
author: "Xiao Xiao"
date: "2021/6/10"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## ScRNAseq data analysis training

Using Dr. Ha's data, here we will show how to do ScRNAseq data analysis with seurat package:

```{r 0}
Sys.setenv(lang="en")   #change error massage into English
setwd("D:/training/scRNAseq")

# 0) Load necessary packages
#install.packages('Seurat')
#install.packages("Matrix", repos="http://R-Forge.R-project.org")
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(openxlsx)
library(ggplot2)
```

## 1) Read in matrix files

Get the barcodes, matrix, and features file from the cell ranger. Then put them together in a folder.

```{r 1}
matrix_dir <- "E:/coding materials/single cell RNA sequencing/210114"
pbmc.data<-Read10X(data.dir = matrix_dir)
```
### 1.1) If you want, you can take a look at the expression data
```{r 1.1}
exprs <- as.matrix(pbmc.data)
```

## 2) Quality control (filter)
### 2.1) Filter cells that has less than 200 gene expression values, and genes that are expressed in less than 3 cells
```{r 2.1}

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```
## 2.2) Filter cells that exhibit aberrantly high gene counts and cells that exhibit extensive mitochondrial contamination (dying cells)

### Mark all the mitochondrial genes
```{r 2.2}
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```
### Draw a violin plot to determine the filtering threshold
``` {r 2.2.1}
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
![image](https://user-images.githubusercontent.com/72723852/121457830-4823d880-c9e3-11eb-98d1-70a9f6cab612.png)

### Apply the filtering threshold
```{r 2.2.2}
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
```
## 2.3) If you want, you can take a look at the expression
```{r 2.3}
exprs2<-as.matrix(pbmc@assays$RNA@data)
```

## 3) Normalize the data 
### 3.1) For each cell, the feature counts are divided by total counts and multiplied by scale.factor. Then the counts are log transformed.
### To take a look at total gene counts before normalization.  count_sum_before<-1:4400;    for (i in 1:4400)  {count_sum_before[i]<-sum(exprs2[,i])};  plot(1:4400,count_sum_before)

```{r 3.1}

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```
### To take a look at total gene counts after normalization.   count_sum_after<-1:4400;  exprs3<-as.matrix(pbmc@assays$RNA@data);     for (i in 1:4400)  {count_sum_after[i]<-sum(exprs3[,i])};  plot(1:4400,count_sum_after)

### 3.2) Before clustering, find the 2000 most variable genes used for clustering
```{r 3.2}
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```
### 3.3) Scale data. This is a standard preprocessing step before PCA. It shifts the expression of each gene, so that the mean expression across cells is 0 - Scales the expression of each gene, so that the variance across cells is 1.
```{r 3.3}
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

```

## 4) Clustering 
### 4.1) Run PCA to reduce the 2000 dimentions to around 10. Here we can only use PCA to reduce dimentions for clustering in seurat.
```{r  4.1}

set.seed(123)
pbmc_2000 <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```
## 4.2) Clustering via k-nearest-neighbor algorithm, using the first 10 PCA dimentions
### In this step, the smaller the k and the larger the resolution, the more clusters we will get.
```{r 4.2}
pbmc_2000_k15 <- FindNeighbors(pbmc_2000, k.param = 15, dims = 1:10)
pbmc_2000_k15_r1 <- FindClusters(pbmc_2000_k15, resolution = 1.0)
```
## 4.3) Run umap for visualization
```{r  4.3}
pbmc_2000_k15_r1_umap <- RunUMAP(pbmc_2000_k15_r1, dims = 1:10)
DimPlot(pbmc_2000_k15_r1_umap, reduction = "umap")
```
![image](https://user-images.githubusercontent.com/72723852/121457881-5f62c600-c9e3-11eb-8f1d-96a34fcac0ed.png)

## 4.4) Decide markers for certain cell types
```{r  4.4}
markers.to.plot <- c("CD3D","CD4","CTLA4","IL7R","CCR7","S100A4","CD14","LYZ","CD19","MS4A1","CD79A","CD79B","BLNK","CD8A","GZMB","CD8B","FCGR3A","MS4A7","NCAM1","KLRB1","KLRC1","KLRD1","KLRF1","GNLY","NKG7","IL3RA","CLEC4C","NRP1","FCER1A","CST3","PPBP")
```
## 4.5) Draw a dot plot to assign clusters
```{r   4.5}
DotPlot(object = pbmc_2000_k15_r1_umap, features = rev(x = markers.to.plot), cols = c("blue", "red"),dot.scale = 8) + RotatedAxis()
```
![image](https://user-images.githubusercontent.com/72723852/121457911-6a1d5b00-c9e3-11eb-9ff9-a78834078ad8.png)

## 4.6) Input your assignment results
```{r  4.6}
pbmc_2000_k15_r1_umap_labelled <- RenameIdents(pbmc_2000_k15_r1_umap, `0` = "Naive CD4+ T cells", `1` = "Naive CD4+ T cells", `2` = "Memory CD4+ T cells", `3` = "Monocytes", `4` = "B cells", `5` = "IL 7R-CD8+ T cells", `6` ="NK cells" , `7` ="NK cells" , `8` = "Memory CD4+ T cells", `9` = "Naive CD8+ T cells", `10` = "Naive CD8+ T cells", `11` = "Naive CD4+ T cells", `12` = "pDC or platelet")
# 4.7) Draw the plot again
DimPlot(pbmc_2000_k15_r1_umap_labelled, label = TRUE)
```
![image](https://user-images.githubusercontent.com/72723852/121457934-71dcff80-c9e3-11eb-9959-f2c27760335a.png)

## 4.8) We can check certain cell marker expression
```{r  4.8}
FeaturePlot(pbmc_2000_k15_r1_umap_labelled, features = c("CD14", "LYZ", "FCGR3A", "MS4A7", "MS4A1"))
```
![image](https://user-images.githubusercontent.com/72723852/121457957-7acdd100-c9e3-11eb-8130-9a49cb99e324.png)

## 4.9) Input the cellular assignment results into the data, for later DEG analysis
```{r  4.9}
pbmc_2000_k15_r1_umap_labelled@meta.data[,6]<-data.frame(pbmc_2000_k15_r1_umap_labelled@active.ident)
```
## 4.10) Get the cell population
### Extract metadata
```{r  4.10}

metadata<-pbmc_2000_k15_r1_umap_labelled@meta.data
```
### Get cell population for control&treat. In this case we know the first 2325 cells are control sample so I just directly write the corresponding code. If in another case the cell number is not 2325, the code will need to be modified accordingly.
```{r  4.10.1}

cell_pop<-data.frame(table(metadata[1:2325,6]),table(metadata[2326:4507,6]))
```
### Delete one extra column and rename the columns
```{r   4.10.2}

cell_pop<-cell_pop[,-3]
colnames(cell_pop)<-c("cell_type","cell_number_control","cell_number_test")
```

### write the cell population into excel file
```{r  4.10.3}

write.xlsx(cell_pop,"cell_pop.xlsx",colnames=TRUE)

```

##  5) DEG analysis
## 5.1) Choose the cell type you wish to analysis. Here we use B cells as an example
```{r 5.1}

B_cells<- subset(pbmc_2000_k15_r1_umap_labelled,subset=seurat_clusters=="B cells")
```
## 5.2) Separate control sample and treatment sample
```{r 5.2}

cellID<-colnames(B_cells@assays$RNA@counts)
sampleID <- as.factor(gsub(".*-","",cellID))
B_cells[["sampleID"]] <- sampleID
B_cells <- SetIdent(object = B_cells, value = "sampleID")
```
## 5.3) Find the DEG using DESeq2 method. Here ident.1 is control sample, ident.2 is treatment sample

```{r 5.3}

deg_B <- FindMarkers(B_cells, ident.1 = 1, ident.2 = 2, logfc.threshold=0.2, min.pct = 0.05, test.use = "DESeq2")
```
## 5.4) Write the DEG results into excel file
```{r 5.4}

write.xlsx(deg_B,"DEG_B.xlsx",row.names=TRUE,colnames=TRUE)
```
## 5.5) Plot heatmap (Here we plot the top 20 DEGs for B cells)
```{r 5.5}

DoHeatmap(B_cells,features = row.names(deg_B)[1:20],group.by = "ident")

```
![image](https://user-images.githubusercontent.com/72723852/121458071-a94bac00-c9e3-11eb-9365-e6832098d532.png)
