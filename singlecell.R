library(Seurat)
library(tidyverse)

## lode the data ##

COL07 <-  ReadMtx(mtx = "separated_samples/COL07/matrix.mtx",
                      features = "separated_samples/COL07/features.tsv",
                      cells = "separated_samples/COL07/barcodes.tsv")

COL12 <-  ReadMtx(mtx = "separated_samples/COL12/matrix.mtx",
                  features = "separated_samples/COL12/features.tsv",
                  cells = "separated_samples/COL12/barcodes.tsv")


COL15 <-  ReadMtx(mtx = "separated_samples/COL15/matrix.mtx",
                  features = "separated_samples/COL15/features.tsv",
                  cells = "separated_samples/COL15/barcodes.tsv")

COL16 <-  ReadMtx(mtx = "separated_samples/COL16/matrix.mtx",
                  features = "separated_samples/COL16/features.tsv",
                  cells = "separated_samples/COL16/barcodes.tsv")


COL17 <-  ReadMtx(mtx = "separated_samples/COL17/matrix.mtx",
                  features = "separated_samples/COL17/features.tsv",
                  cells = "separated_samples/COL17/barcodes.tsv")

COL18 <-  ReadMtx(mtx = "separated_samples/COL18/matrix.mtx",
                  features = "separated_samples/COL18/features.tsv",
                  cells = "separated_samples/COL18/barcodes.tsv")

#Create seurat Object

COL07 <- CreateSeuratObject(counts = COL07 ,project = "COL07",min.cells = 3,min.features = 200)

COL07 <- PercentageFeatureSet(COL07, pattern = "^MT-", col.name = "percent.mt")

COL12 <- CreateSeuratObject(counts = COL12 ,project = "COL12",min.cells = 3,min.features = 200)

COL12 <- PercentageFeatureSet(COL12, pattern = "^MT-", col.name = "percent.mt")

COL15 <- CreateSeuratObject(counts = COL15 ,project = "COL15",min.cells = 3,min.features = 200)

COL15 <- PercentageFeatureSet(COL15, pattern = "^MT-", col.name = "percent.mt")

COL16 <- CreateSeuratObject(counts = COL16 ,project = "COL16",min.cells = 3,min.features = 200)

COL16 <- PercentageFeatureSet(COL16, pattern = "^MT-", col.name = "percent.mt")

COL17 <- CreateSeuratObject(counts = COL17 ,project = "COL17",min.cells = 3,min.features = 200)

COL17 <- PercentageFeatureSet(COL17, pattern = "^MT-", col.name = "percent.mt")

COL18 <- CreateSeuratObject(counts = COL18 ,project = "COL18",min.cells = 3,min.features = 200)

COL18 <- PercentageFeatureSet(COL18, pattern = "^MT-", col.name = "percent.mt")



#filteration

COL07 <- subset(COL07, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10)
COL12 <- subset(COL12, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10 )

COL15 <- subset(COL15, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10 )
COL16 <- subset(COL16, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10 )

COL17 <- subset(COL17, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10 )
COL18 <- subset(COL18, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10 )


# Combine the two datasets into a list

datasets_list <- list(COL07 = COL07, COL12 = COL12, COL15 = COL15, COL16 = COL16, COL17 = COL17, COL18 = COL18)

rm(COL07, COL12, COL15, COL16, COL17, COL18)

# normalize and identify variable features for eatch data set independentaly

datasets_list <- lapply( X = datasets_list , FUN = function(X) {
  X = NormalizeData(X)
  X = FindVariableFeatures( X , selection.method = "vst" , nfeatures = 2000)
})


#select integration features

features <- SelectIntegrationFeatures(object.list = datasets_list)

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = datasets_list, anchor.features = features)

# Find integration anchors and integrate the data
#integration_anchors <- FindIntegrationAnchors(object.list = list(early_tumor_1, late_tumor_1), dims = 1:30)
integrated_data <- IntegrateData(anchorset = integration_anchors)


COL07@meta.data$percent.mt

unique(GSE249002@meta.data$percent.mt)

COL07@meta.data

#filtering low quality cells

VlnPlot(GSE249002, features = c( "nCount_RNA","nFeature_RNA" ,"percent.mt"),ncol = 3)
GSE249002 <- subset(GSE249002, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10)

FeatureScatter(GSE249002, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

view(GSE249002@meta.data)
## normalize theta

GSE249002 <- NormalizeData(GSE249002)

str(GSE249002)
##identify the highly variable features

GSE249002 <- FindVariableFeatures(GSE249002, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(GSE249002), 10)

##5.scaling the data
all.genes <- rownames(GSE249002)

GSE249002 <- ScaleData(GSE249002, features = all.genes )

## perfome linear dimensionally reduction

GSE249002 <- RunPCA(GSE249002, features = VariableFeatures(GSE249002))
DimHeatmap(GSE249002,dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(GSE249002)

## clustring

GSE249002 <- FindNeighbors(object = GSE249002 , reduction = "pca", dims = 4:15)

GSE249002 <- FindClusters(object = GSE249002,resolution =  c(0.3))
view(GSE249002@meta.data)
colnames(GSE249002@meta.data)
unique(GSE249002@meta.data$RNA_snn_res.0.3)

DimPlot(GSE249002, group.by = "RNA_snn_res.0.3", label = TRUE)

## Setting identity of the clusters
Idents(GSE249002)
Idents(GSE249002) <- "RNA_snn_res.0.1"

## non linear dimentionally reduction

GSE249002 <- RunUMAP(object = GSE249002, reduction = "pca",dims = 1:15 )

DimPlot(GSE249002, reduction ="umap", label = TRUE)

saveRDS(GSE249002, file = "GSE249002.rds")


#install.packages("DoubletFinder")

## annotations by using SingleR

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(SingleCellExperiment)

ref_mouse <- celldex::MouseRNAseqData()

View(as.data.frame(colData(ref_mouse)))

sce <- as.SingleCellExperiment(GSE249002)

annotations_mouse <- SingleR(test = sce, ref = ref_mouse, labels = ref_mouse$label.main)

GSE249002$SingleR_RNAseq_labels <- annotations_mouse$labels

DimPlot(GSE249002, reduction = "umap", group.by = "SingleR_RNAseq_labels", label = TRUE)

DimPlot(GSE249002, reduction = "umap", label = TRUE)

view(GSE249002@meta.data)



#Identify all marker genes
markers <- FindAllMarkers(GSE249002, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'DESeq2',slot = 'counts')


top_markers_cluster0 <- markers %>% filter(cluster == 0) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster0)


FeaturePlot(GSE249002, features = "Trim47", min.cutoff = "q5")


GSE249002 <- RenameIdents(GSE249002, '0' = "B cells")
Idents(GSE249002)


top_markers_cluster1 <- markers %>% filter(cluster == 1) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster1)

GSE249002 <- RenameIdents(GSE249002, `1` = " T cell cluster")


top_markers_cluster2 <- markers %>% filter(cluster == 2) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster2)
unique(markers$cluste)

GSE249002 <- RenameIdents(GSE249002, `2` = "B cells")

top_markers_cluster3 <- markers %>% filter(cluster == 3) %>% top_n(n = 20, wt = avg_log2FC)
print(top_markers_cluster3)

GSE249002 <- RenameIdents(GSE249002, `3` = "Granulocytes")

top_markers_cluster4 <- markers %>% filter(cluster == 4) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster4)

GSE249002 <- RenameIdents(GSE249002, `4` = "T follicular helper (Tfh) cells")

top_markers_cluster5 <- markers %>% filter(cluster == 5) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster5)

GSE249002 <- RenameIdents(GSE249002, `5` = "Monocytes")

top_markers_cluster6 <- markers %>% filter(cluster == 6) %>% top_n(n = 20, wt = avg_log2FC)
print(top_markers_cluster6)

GSE249002 <- RenameIdents(GSE249002, `6` = "Effector CD8+ T Cells")

top_markers_cluster7 <- markers %>% filter(cluster == 7) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster7)

GSE249002 <- RenameIdents(GSE249002, `7` = "Activated or memory T ")

VlnPlot(GSE249002, features = c("Tmem123", "Lsp1", "Epsti11", "Crip12", "Bst21", "Tbc1d4"))


top_markers_cluster8 <- markers %>% filter(cluster == 8) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster8)

GSE249002 <- RenameIdents(GSE249002, `8` = "Proliferating Cells")

top_markers_cluster9 <- markers %>% filter(cluster == 9) %>% top_n(n = 10, wt = avg_log2FC)
print(top_markers_cluster9)

GSE249002 <- RenameIdents(GSE249002, `9` = "Epithelial Cells")

table(Idents(GSE249002), GSE249002$SingleR_RNAseq_labels)



cluster_4_cells <- WhichCells(GSE249002, idents = 4)
genes_in_cluster_4 <- rownames(GSE249002)[which(GSE249002@meta.data[cluster_4_cells, "seurat_clusters"] == 4)]

cluster_1_cells <- WhichCells(GSE249002, idents = 1)
genes_in_cluster_1 <- rownames(GSE249002)[which(GSE249002@meta.data[cluster_1_cells, "seurat_clusters"] == 4)]
print(genes_in_cluster_4)


DimPlot(GSE249002, reduction ="umap", label = TRUE)

VlnPlot(GSE249002, features = c("Cd8a", "Cd8b1"))

write.csv( markers, file = "deg_GSE249002.csv", row.names = TRUE)

saveRDS(GSE249002 ,"GSE249002.rds")


GSE249002$seurat_clusters <- Idents(GSE249002)

Tcell_EffectorTcell <- FindMarkers(GSE249002, ident.1 = " T cell cluster", ident.2 = "Effector CD8+ T Cells",
                                   logfc.threshold = 0, test.use = "wilcox")

write.csv(Tcell_EffectorTcell, file = "deg_Tcell_EffectorTcell.csv", row.names = TRUE)
