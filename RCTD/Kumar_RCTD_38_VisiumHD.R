library(Seurat)
library(spacexr)
library(SingleCellExperiment)
library(Matrix)

setwd("PATH")
# ----------------------------
# Step 1: Load scRNA Seurat object
# ----------------------------
load("snRNAseq_subset.rda")

counts <- snRNAseq_subset@assays$RNA@layers$counts
rownames(counts) <- rownames(snRNAseq_subset)
colnames(counts) <- colnames(snRNAseq_subset)


cell_types <- snRNAseq_subset$cell_type
names(cell_types) <- colnames(snRNAseq_subset) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- snRNAseq_subset$nCount_RNA; names(nUMI) <- colnames(snRNAseq_subset) # create nUMI named list

reference <- Reference(counts, cell_types, nUMI)
# ----------------------------
# Step 2: Load VisiumHD Seurat object
# ----------------------------
load("Patient_38_ER_RCTD.rda")

# Get tissue coordinates using Seurat v5 function
coords <- GetTissueCoordinates(Patient_38_RCTD)

# Extract counts matrix
counts <- as.matrix(Patient_38_RCTD@assays$Spatial$counts)

# Sanity check
stopifnot(all(rownames(coords) == colnames(counts)))

# Create SpatialRNA object for RCTD
spatial <- SpatialRNA(coords = coords[, c("x", "y")], counts = counts)

# ----------------------------
# Step 3: Run RCTD
# ----------------------------
rctd <- create.RCTD(spatialRNA = spatial, reference = reference, max_cores = 4)

# Run the model (doublet mode is recommended for complex tissues)
rctd <- run.RCTD(rctd, doublet_mode = 'full')
Patient_38_RCTD <- subset(Patient_38_RCTD, cells = rownames(rctd@results$weights))
# ----------------------------
# Step 4: Add results to Seurat object
# ----------------------------
results <- rctd@results$weights  # cell type proportions per spot
Kumar_celltype <- as.data.frame(results)

# Add as metadata to Seurat object
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Kumar_celltype)
Kumar_dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
Kumar_dominant_df <- data.frame(Kumar_dominant_celltype)
rownames(Kumar_dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Kumar_dominant_df)

pdf("Patient_38_Bhupinder_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="Kumar_dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

save(Patient_38_RCTD, file="Patient_38_RCTD_snRNAseq.rda")


# Adding metadata according to Buphinder and Kumar weights:

setwd("PATH")

############ Patient 38 ###################################

load("Patient_38_RCTD_snRNAseq.rda")


results <- data.frame(
  B_cells = Patient_38_RCTD$`B cells`,
  Endothelial_cells = Patient_38_RCTD$`Endothelial cells`,
  Epithelial_cells = Patient_38_RCTD$`Epithelial cells`,
  Malignant = Patient_38_RCTD$malignant,
  Myeloid_Buphinder = Patient_38_RCTD$Myeloid,
  Plasmablasts = Patient_38_RCTD$plasmablasts,
  Stromal = Patient_38_RCTD$Stromal,
  T_cells_Buphinder = Patient_38_RCTD$`T cells`,
  Lymphatic = Patient_38_RCTD$Lymphatic,
  Adipocytes = Patient_38_RCTD$Adipocytes,
  Vascular = Patient_38_RCTD$Vascular,
  Myeloid_Kumar = Patient_38_RCTD$Myel_Kumar,
  Fibroblast = Patient_38_RCTD$Fibroblast,
  Perivascular = Patient_38_RCTD$Perivascular,
  Mast = Patient_38_RCTD$Mast,
  T_cells_Kumar = Patient_38_RCTD$T_cells,
  Basal = Patient_38_RCTD$Basal,
  LumSec = Patient_38_RCTD$LumSec,
  LumHR = Patient_38_RCTD$LumHR
)

# Add as metadata to Seurat object
Final_dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
Final_dominant_df <- data.frame(Final_dominant_celltype)
rownames(Final_dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Final_dominant_df)


pdf("Patient_38_final_dominat_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="Final_dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

Idents(Patient_38_RCTD) <- Patient_38_RCTD$Final_dominant_celltype

Final_dominant_celltype_ptx38 <- FindAllMarkers(Patient_38_RCTD, only.pos = TRUE)

write.table(Final_dominant_celltype_ptx38, file = "Final_dominant_celltype_ptx38_markers.txt", sep = "\t")

save(Patient_38_RCTD, file="Patient_38_RCTD_snRNAseq.rda")


Final_dominant_celltype_ptx38 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_final_dominant_celltype

pdf("Patient_38_top10markers_dotplot.pdf", height = 7, width = 30)
DotPlot(Patient_38_RCTD, features = unique(top10_final_dominant_celltype$gene)) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
dev.off()

