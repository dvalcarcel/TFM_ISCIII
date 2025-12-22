library(Seurat)
library(spacexr)
library(SingleCellExperiment)
library(Matrix)

# ----------------------------
# Step 1: Load scRNA Seurat object
# ----------------------------
load("PATH/Bhupinder_2021.processed.seurat.RData")  

counts <- processed.seurat@assays$RNA@counts
rownames(counts) <- rownames(processed.seurat)
colnames(counts) <- colnames(processed.seurat)


cell_types <- processed.seurat$celltype
names(cell_types) <- colnames(processed.seurat) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- processed.seurat$nCount_RNA; names(nUMI) <- colnames(processed.seurat) # create nUMI named list

reference <- Reference(counts, cell_types, nUMI)
# ----------------------------
# Step 2: Load VisiumHD Seurat object
# ----------------------------
load("PATH/Patient_38_16um.rda")  

# Get tissue coordinates using Seurat v5 function
coords <- GetTissueCoordinates(Patient_38)

# Extract counts matrix
counts <- as.matrix(Patient_38@assays$Spatial$counts)

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
Patient_38_RCTD <- subset(Patient_38, cells = rownames(rctd@results$weights))
# ----------------------------
# Step 4: Add results to Seurat object
# ----------------------------
results <- rctd@results$weights  # cell type proportions per spot
Bhupinder_celltype <- as.data.frame(results)

# Add as metadata to Seurat object
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = Bhupinder_celltype)
dominant_celltype <- colnames(results)[apply(results, 1, which.max)]
dominant_df <- data.frame(dominant_celltype)
rownames(dominant_df) <- rownames(results)
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = dominant_df)

setwd("PATH/RCTD")

pdf("Patient_38_Bhupinder_celltype.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="dominant_celltype", image.alpha=0, pt.size.factor=3.2)
dev.off()

save(Patient_38_RCTD, file="Patient_38_RCTD.rda")
