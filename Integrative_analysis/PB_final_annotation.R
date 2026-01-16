library(Seurat)
library(dplyr)
library(ggplot2)

setwd("PATH")


# Loading all patients

load("../Patient_4_26_RCTD_snRNAseq.rda")
load("../Patient_06_RCTD_snRNAseq.rda")
load("../Patient_38_RCTD_snRNAseq.rda")
load("../Patient_41_67_RCTD_snRNAseq.rda")
load("../Patient_61_RCTD_snRNAseq.rda")

Patient_4_26_RCTD <- Patient_4_26_RCTD[, !is.na(Patient_4_26_RCTD$final_annotation)]
Patient_06_RCTD <- Patient_06_RCTD[, !is.na(Patient_06_RCTD$final_annotation)]
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$final_annotation)]
Patient_41_67_RCTD <- Patient_41_67_RCTD[, !is.na(Patient_41_67_RCTD$final_annotation)]
Patient_61_RCTD <- Patient_61_RCTD[, !is.na(Patient_61_RCTD$final_annotation)]

# We make a data.frame for each patient and merge it
#  Patient 6 block a1, 38 block A2, 4_26 block B1, 61 block B2, 41_67 block C1


df_06 <- data.frame(
  spot = paste0("A1_", names(Patient_06_RCTD$final_annotation)),
  final_annotation = Patient_06_RCTD$final_annotation
)

df_38 <- data.frame(
  spot = paste0("A2_", names(Patient_38_RCTD$final_annotation)),
  final_annotation = Patient_38_RCTD$final_annotation
)

df_4_26 <- data.frame(
  spot = paste0("B1_", names(Patient_4_26_RCTD$final_annotation)),
  final_annotation = Patient_4_26_RCTD$final_annotation
)

df_61 <- data.frame(
  spot = paste0("B2_", names(Patient_61_RCTD$final_annotation)),
  final_annotation = Patient_61_RCTD$final_annotation
)

df_41_67 <- data.frame(
  spot = paste0("C1_", names(Patient_41_67_RCTD$final_annotation)),
  final_annotation = Patient_41_67_RCTD$final_annotation
)

final_annotation_df <- rbind(df_4_26, df_06, df_38, df_41_67, df_61)


# Assuring that all spots are as rownames and deleting spot column
rownames(final_annotation_df) <- final_annotation_df$spot

final_annotation_df$spot <- NULL



# Loading integrated object and adding the metadata 

load("ST_ER_RCTD_integrated.rda")

ST_integrated <- AddMetaData(ST_integrated, metadata = final_annotation_df)
ST_integrated$final_annotation_Denosumab <-paste(ST_integrated$final_annotation, ST_integrated$Denosumab, sep = "_")

Idents(ST_integrated) <- ST_integrated$final_annotation_Denosumab

save(ST_integrated, file = "ST_ER_RCTD_integrated.rda")


# UMAPs

# Defining the same color for each annotation as SpatialDimPlots

cols <- c(
  ductal_invasion        = "#4A0000",  
  tumor                  = "#7F0000",  
  tumor_mixed            = "#F70000",
  "in situ carcinoma"    = "#2A0000",
  normal_breast_tissue   = "#1A9850", 
  normal_breast_mixed    = "#66C2A5",  
  Stromal                = "#BCA9FC",  
  stroma                 = "#BCA9FC",
  artery                 = "#4126AC",
  Endothelial_cells      = "#2166AC",  
  Lymphatic              = "#4393C3", 
  Perivascular           = "#92C5DE",  
  B_cells                = "#F46D43",  
  T_cells                = "#FDAE61",  
  Plasmablasts           = "#1CECFC",  
  Myeloid                = "#ED91E5",  
  Mast                   = "#FCE99D",  
  Adipocytes             = "#FFFF33"   
)


p1 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Patient") +
  ggtitle("A) Patient") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16)) 
p2 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, group.by = "Denosumab") +
  ggtitle("B) Denosumab") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16)) 
p3 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, cols = cols, group.by = "Pathologist_annotation") +
  ggtitle("C) Pathologist annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16))
p4 <- DimPlot(object = ST_integrated, reduction = "umap", raster=FALSE, cols = cols, group.by = "final_annotation") +
  ggtitle("D) Final annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16))

png("UMAPs_ST_integrated.png", width = 12000, height = 10000, res = 650)
p1+p2+p3+p4
dev.off()

# Counts matrix with the final annotation for each patient

agg_expr <- AggregateExpression(
  object = ST_integrated,
  group.by = c("final_annotation_Denosumab", "Patient"),
  assay = "Spatial",
  slot = "counts"  
)

pb_matrix <- as.matrix(agg_expr$Spatial)
write.table(pb_matrix, file = "ST_integrated_pb_matrix_final_annotation.txt", sep = '\t')

#FindAllmarkers of all idents

ST_integrated.markers <- FindAllMarkers(ST_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(ST_integrated.markers, "ST_integrated_markers_final_annotation.txt", sep = '\t')

#Filtering top5 markers

ST_integrated.markers <- read.table("ST_integrated_markers_final_annotation.txt", header = TRUE, sep = '\t')

ST_integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5_final_dominant_celltype

# To get unique clusters
clusters <- unique(ST_integrated.markers$cluster)

# To split in Pre and Post
clusters_pre  <- clusters[grepl("_Pre$", clusters)]
clusters_post <- clusters[grepl("_Post$", clusters)]

# To delete suffix "_Pre" and "_Post" 
base_pre  <- gsub("_Pre$", "", clusters_pre)
base_post <- gsub("_Post$", "", clusters_post)

# Alphabetical sorting
base_order <- sort(unique(c(base_pre, base_post)))

# Final cluster order
orden_clusters <- unlist(lapply(base_order, function(x) {
  c(paste0(x, "_Pre"), paste0(x, "_Post"))
}))

orden_clusters

# Final order within the idents
ST_integrated$cluster <- factor(Idents(ST_integrated), 
                                levels = orden_clusters)
Idents(ST_integrated) <- factor(Idents(ST_integrated), levels = orden_clusters)


# Saving dotplot markers to PDF
pdf("Final_annotation_top5markers_dotplot.pdf", height = 7, width = 40)
DotPlot(ST_integrated, features = unique(top5_final_dominant_celltype$gene)) +
  ggtitle("Top 5 markers dotplot by cell identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold")) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
dev.off()


# Clean dotplot with Top5 markers of Endothelial, Tumor, T, B and Myeloid cells

idents_to_filter <- c("Endothelial_cells_Post", "T_cells_Post", "tumor_Post", "B_cells_Post", "Myeloid_Post",
                     "tumor_Pre", "Myeloid_Pre", "Endothelial_cells_Pre", "T_cells_Pre", "B_cells_Pre",
                     "Stromal_Pre", "Stromal_Post")

ST_integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(cluster %in% idents_to_filter) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5_final_dominant_celltype_filtered

png("Final_annotation_top5markers_dotplot_filtered.png", width = 7200, height = 2400, res = 450)
DotPlot(ST_integrated, features = unique(top5_final_dominant_celltype_filtered$gene), idents = idents_to_filter) +
  ggtitle("F) Top 5 markers dotplot by cell identity (filtered)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
dev.off()


