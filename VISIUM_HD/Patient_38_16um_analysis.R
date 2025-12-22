### Visium HD analysis from: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette ###

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)

local_dir <- "PATH/square_016um/"
Patient_38 <- Load10X_Spatial(data.dir = local_dir,filename = "filtered_feature_bc_matrix.h5",assay = "Spatial")
# Setting default assay changes between 8um and 16um binning
Assays(Patient_38)
# DefaultAssay(Patient_38) <- "Spatial.008um"

## Add metadata

#Patient_38_metadata <- read.csv('metadata_blockB2.csv', header = TRUE, sep = ',')
#Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Block, col.name = "Block")
#Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Patient, col.name = "Patient")
#Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_metadata$Denosumab, col.name = "Denosumab")


pdf("Patient_38_nCount_16um.pdf", height = 7, width = 14)
vln.plot <- VlnPlot(Patient_38, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(Patient_38, features = "nCount_Spatial",pt.size.factor = 3.2) + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot
dev.off()

# Patient_38 <- subset(Patient_38, subset = nCount_Spatial > 200)

## Normalize data
# normalize both 8um and 16um bins

Patient_38 <- NormalizeData(Patient_38)
Patient_38 <- FindVariableFeatures(Patient_38)
Patient_38 <- ScaleData(Patient_38)
Patient_38 <- RunPCA(Patient_38)
Patient_38 <- RunUMAP(Patient_38, dims = 1:50)

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


Patient_38 <- FindNeighbors(Patient_38, dims = 1:50)
Patient_38 <- FindClusters(Patient_38, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

pdf("clustree_Patient_38_16um.pdf",width = 10, height = 14)
clustree(Patient_38)
dev.off()

Idents(Patient_38) <- Patient_38$Spatial_snn_res.0.7

pdf("Patient_38_clusters_res_0.7.pdf")
SpatialDimPlot(Patient_38, image.alpha = 0,pt.size.factor=3.2)
dev.off()

Sign_names <- read.table("Sign_names.txt", sep = '\t')
Sign_names <- Sign_names$V1

pdf("Dotplot_signatures_clusters.pdf", height = 7, width = 14)
DotPlot(Patient_38, features = Sign_names)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_color_gradientn(colors = c("blue", "grey", "red"))
dev.off()

Patient_38_0.7_markers <- FindAllMarkers(Patient_38, only.pos = TRUE)

Patient_38_0.7_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

Idents(Patient_38) <- factor(Idents(Patient_38), levels = as.character(0:15))


pdf("Dot_plot_top_markers.pdf", height = 7, width = 20)
DotPlot(Patient_38, features = unique(top10$gene)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8)) +
  scale_color_gradientn(colors = c("blue", "grey", "red"))
dev.off()

library(scCATCH)
Patient_38_catch <- createscCATCH(data = Patient_38[["Spatial"]]@layers$counts, cluster = as.character(Patient_38$Spatial_snn_res.0.7))
Patient_38_catch <- findmarkergene(object = Patient_38_catch, species = "Human", marker = cellmatch,cancer = TRUE,tissue = "Breast", use_method = "1")
Patient_38_catch <- findcelltype(HEALTHY_1_catch)

library(VISION) 

Signature_scores <- as.data.frame(Patient_38_vis@SigScores)
all(rownames(Signature_scores) %in% colnames(Patient_38))

Patient_38 <- AddMetaData(Patient_38, metadata = Signature_scores)

# Define output directory for PDFs
output_dir <- "Signature_plots/"
dir.create(output_dir, showWarnings = FALSE)

# Extract signature names from metadata
signature_names <- colnames(Signature_scores)

# Generate and save SpatialFeaturePlots for each signature
for (signature in signature_names) {
  p <- SpatialFeaturePlot(Patient_38, features = signature,image.alpha=0,pt.size.factor=3.2) + ggtitle(signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, ".pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}

print("All spatial feature plots have been saved.")

library(loupeR)

create_loupe_from_seurat(Patient_38, output_name = "Patient_38_16um.cloupe")
Patient_38_metadata <- read.csv("PATH/metadata_blockA2.csv")
Patient_38_16um_bc <- as.data.frame(colnames(Patient_38))
Patient_38_16um_bc$`colnames(Patient_38)` <- gsub("016um", "008um", Patient_38_16um_bc$`colnames(Patient_38)`)
Patient_38_16um_metadata <- Patient_38_metadata[Patient_38_metadata$Barcode %in% Patient_38_16um_bc$`colnames(Patient_38)`, ]

Patient_38_16um_metadata$Barcode <- gsub("008um", "016um", Patient_38_16um_metadata$Barcode)

## Add metadata

Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Block, col.name = "Block")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Patient, col.name = "Patient")
Patient_38 <- AddMetaData(object = Patient_38, metadata = Patient_38_16um_metadata$Denosumab, col.name = "Denosumab")

save(Patient_38, file = "Patient_38_16um.rda")

