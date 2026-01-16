library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(clustree)

# This script was run in the CNIO's cluster

setwd("PATH")
options(future.globals.maxSize= 5*1024*1024^2)

load("../Patient_06_ER_RCTD.rda")
load("../Patient_38_ER_RCTD.rda")
load("../Patient_4_26_ER_RCTD.rda")
load("../Patient_61_ER_RCTD.rda")
load("../Patient_41_67_ER_RCTD.rda")

ST_merged  <- merge(x = Patient_06_RCTD, y = c(Patient_38_RCTD, Patient_4_26_RCTD, Patient_61_RCTD, Patient_41_67_RCTD), add.cell.ids = c("A1","A2","B1","B2","C1"), merge.data = TRUE)

ST_merged <- JoinLayers(ST_merged,  assay = "Spatial")
ST_merged <- subset(ST_merged, subset = nCount_Spatial > 50)

ST_merged <- NormalizeData(ST_merged)
ST_merged <- FindVariableFeatures(ST_merged)
ST_merged <- ScaleData(ST_merged)
ST_merged <- RunPCA(ST_merged)

pdf('elbow_ST_merged.pdf')
ElbowPlot(ST_merged,ndims=50)
dev.off()

ST_merged <- FindNeighbors(ST_merged, dims = 1:30)
ST_merged <- FindClusters(ST_merged,resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf('clustree_ST_merged.pdf',width = 14, height = 14)
clustree(ST_merged)
dev.off()

ST_merged <- RunUMAP(ST_merged,dims = 1:30)

Idents(ST_merged) <- ST_merged$Spatial_snn_res.0.1

pdf('UMAP_ST_merged.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, label = TRUE)
dev.off()

pdf('UMAP_ST_merged_patient.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, group.by = "Patient")
dev.off()

pdf('UMAP_ST_block.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, group.by = "Block")
dev.off()

pdf('UMAP_ST_denosumab.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, group.by = "Denosumab")
dev.off()

pdf('UMAP_ST_Pathologist_annotation.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, group.by = "Pathologist_annotation")
dev.off()

pdf('UMAP_ST_RCTD_Bhupinder_annotation.pdf')
DimPlot(object = ST_merged, reduction = "umap", raster=FALSE, group.by = "dominant_celltype")
dev.off()

Spots_to_remove <- read.table("../../../Integrated_analysis/Spots_to_remove.txt", sep ='\t')
Spots_to_remove <- Spots_to_remove$V1
ST_merged <- subset(ST_merged, cells= Spots_to_remove, invert=TRUE)

save(ST_merged, file="ST_ER_RCTD_merged.rda")

library(harmony)

## Harmony integration:

ST_integrated <- RunHarmony(ST_merged, "Block",kmeans_init_nstart=20, kmeans_init_iter_max=100)

ST_integrated <- RunUMAP(ST_integrated, reduction="harmony", dims=1:30)

pdf('UMAP_ST_integrated_patient.pdf')
DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Patient")
dev.off()

pdf('UMAP_ST_integrated_block.pdf')
DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Block")
dev.off()

pdf('UMAP_ST_integrated_denosumab.pdf')
DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Denosumab")
dev.off()

pdf('UMAP_ST_integrated_Pathologist_annotation.pdf')
DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Pathologist_annotation")
dev.off()

pdf('UMAP_ST_integrated_RCTD_Bhupinder_annotation.pdf')
DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "dominant_celltype")
dev.off()

ST_integrated <- FindNeighbors(ST_integrated, dims = 1:30)
ST_integrated <- FindClusters(ST_integrated,resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf('clustree_ST_integrated.pdf',width = 14, height = 14)
clustree(ST_integrated)
dev.off()

ST_integrated.markers <- FindAllMarkers(ST_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(ST_integrated.markers, "ST_integrated_markers.txt", sep = '\t')


save(ST_integrated, file="ST_ER_RCTD_integrated.rda")
