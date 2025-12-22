### Visium HD analysis from: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette ###

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)

local_dir <- "PATH"
Patient_38 <- Load10X_Spatial(data.dir = local_dir,filename = "filtered_feature_bc_matrix.h5",assay = "Spatial")
# Setting default assay changes between 8um and 16um binning
Assays(Patient_38)
# DefaultAssay(Patient_38) <- "Spatial.008um"

## Add metadata

Patient_38_metadata <- read.csv('metadata_blockA2.csv', header = TRUE, sep = ',')
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Block, col.name = "Block")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Patient, col.name = "Patient")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Denosumab, col.name = "Denosumab")


pdf("Patient_38_nCount.pdf", height = 7, width = 14)
vln.plot <- VlnPlot(Patient_38, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(Patient_38, features = "nCount_Spatial",pt.size.factor = 3.2) + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot
dev.off()

## Normalize data
# normalize both 8um and 16um bins

Patient_38 <- NormalizeData(Patient_38)

## Visualize gene expression

pdf("Patient_38_EPCAM_PTPRC.pdf", height = 7, width = 14)
p1 <- SpatialFeaturePlot(Patient_38, features = "EPCAM", image.alpha = 0,pt.size.factor = 3.2) + ggtitle("EPCAM")
p2 <- SpatialFeaturePlot(Patient_38, features = "PTPRC",image.alpha = 0,pt.size.factor = 3.2) + ggtitle("PTPRC")
p1 | p2
dev.off()

pdf("Patient_38_ESR1_IGKC.pdf", height = 7, width = 14)
p1 <- SpatialFeaturePlot(Patient_38, features = "ESR1", image.alpha = 0,pt.size.factor = 3.2) + ggtitle("ESR1")
p2 <- SpatialFeaturePlot(Patient_38, features = "IGKC",image.alpha = 0,pt.size.factor = 3.2) + ggtitle("IGKC")
p1 | p2
dev.off()

## Unsupervised clustering

# note that data is already normalized
Patient_38 <- FindVariableFeatures(Patient_38)
Patient_38 <- ScaleData(Patient_38)
# we select 50,0000 cells and create a new 'sketch' assay
Patient_38 <- SketchData(
  object = Patient_38,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(Patient_38) <- "sketch"

# perform clustering workflow
Patient_38 <- FindVariableFeatures(Patient_38)
Patient_38 <- ScaleData(Patient_38)
Patient_38 <- RunPCA(Patient_38, assay = "sketch", reduction.name = "pca.sketch")
Patient_38 <- FindNeighbors(Patient_38,assay = "sketch", reduction = "pca.sketch", dims = 1:50)
Patient_38 <- FindClusters(Patient_38, cluster.name = "seurat_cluster.sketched", resolution = c(0.5))

# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# Number of nodes: 50000
# Number of edges: 1741909

# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.8541
#    Number of communities: 113
#    Elapsed time: 30 seconds
#    96 singletons identified. 17 final clusters.

Patient_38 <- RunUMAP(Patient_38, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

## Projecting sketch analysis into the entire dataset:

Patient_38 <- ProjectData(
  object = Patient_38,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(Patient_38) <- "sketch"
Idents(Patient_38) <- "seurat_cluster.sketched"
p1 <- DimPlot(Patient_38, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(Patient_38) <- "Spatial"
Idents(Patient_38) <- "seurat_cluster.projected"
p2 <- DimPlot(Patient_38, reduction = "full.umap.sketch", label = F,raster = FALSE) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

pdf("Patient_38_clustering.pdf", height = 10, width = 14)
p1 | p2
dev.off()

pdf("Patient_38_spatial_clustering.pdf")
SpatialDimPlot(Patient_38, label = T, repel = T, label.size = 4,image.alpha = 0,pt.size.factor = 3.2)
dev.off()

pdf("Patient_38_spatial_block.pdf")
SpatialDimPlot(Patient_38,label.size = 4,image.alpha = 0,pt.size.factor = 3.2,group.by = "Block")
dev.off()

pdf("Patient_38_spatial_patient.pdf")
SpatialDimPlot(Patient_38,label.size = 4,image.alpha = 0,pt.size.factor = 3.2,group.by = "Patient")
dev.off()

pdf("Patient_38_spatial_denosumab.pdf")
SpatialDimPlot(Patient_38,label.size = 4,image.alpha = 0,pt.size.factor = 3.2,group.by = "Denosumab")
dev.off()

# Create downsampled object to make visualization either
DefaultAssay(Patient_38) <- "Spatial"
Idents(Patient_38) <- "seurat_cluster.projected"
Patient_38_subset <- subset(Patient_38, cells = Cells(Patient_38[["Spatial"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(Patient_38_subset) <- "Spatial"
Idents(Patient_38_subset) <- "seurat_cluster.projected"
Patient_38_subset <- BuildClusterTree(Patient_38_subset, assay = "Spatial", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(Patient_38_subset, assay = "Spatial", only.pos = TRUE)
write.table(markers, file="Patient_38_markers_0.5.txt", sep = '\t')
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

Patient_38_subset <- ScaleData(Patient_38_subset, assay = "Spatial", features = top5$gene)
p <- DoHeatmap(Patient_38_subset, assay = "Spatial", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()

pdf("Patient_38_HM_top5.pdf", height = 14, width = 10)
p
dev.off()

### Identifying spatially-defined tissue domains

library(SeuratWrappers)
library(Banksy)

Patient_38 <- RunBanksy(Patient_38,
                        lambda = 0.8, verbose = TRUE,
                        assay = "Spatial", slot = "data", features = "variable",
                        k_geom = 50
)

DefaultAssay(Patient_38) <- "BANKSY"
Patient_38 <- RunPCA(Patient_38, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(Patient_38), npcs = 30)
Patient_38 <- FindNeighbors(Patient_38, reduction = "pca.banksy", dims = 1:30)
Patient_38 <- FindClusters(Patient_38, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(Patient_38) <- "banksy_cluster"

pdf('Patient_38_Banksy_clusters.pdf')
p <- SpatialDimPlot(Patient_38, group.by = "banksy_cluster", label = T, repel = T, label.size = 4,image.alpha = 0,pt.size.factor = 3.2)
p
dev.off()

markers <- FindAllMarkers(Patient_38, assay = "BANKSY", only.pos = TRUE)
write.table(markers, file="Patient_38_BANKSY_markers_0.5.txt", sep = '\t')

banksy_cells <- CellsByIdentities(Patient_38)
pdf('Patient_38_Banksy_cells.pdf')
p <- SpatialDimPlot(Patient_38, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T,image.alpha = 0,pt.size.factor = 3.2) + NoLegend()
p
dev.off()

save(Patient_38, file = "Patient_38.rda")
