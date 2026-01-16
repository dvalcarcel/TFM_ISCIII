library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(jsonlite)

#This script was run on the CNIO's cluster

setwd("PATH")

load("Patient_38_RCTD_snRNAseq.rda")
# Remove spots where Denosumab is NA
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$Denosumab)]
Patient_38_RCTD <- Patient_38_RCTD[, !is.na(Patient_38_RCTD$final_annotation)]
Patient_38_RCTD$final_annotation_Denosumab <- paste(Patient_38_RCTD$final_annotation,Patient_38_RCTD$Denosumab,sep="_")
Idents(Patient_38_RCTD) <- Patient_38_RCTD$final_annotation_Denosumab

# Extract expression data and metadata
data.input <- GetAssayData(Patient_38_RCTD, assay = "Spatial", slot = "data")
meta <- data.frame(labels = Idents(Patient_38_RCTD), row.names = names(Idents(Patient_38_RCTD)))
spatial.locs = GetTissueCoordinates(Patient_38_RCTD, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs$cell <- NULL

scalefactors = fromJSON(txt = file.path("PATH", 'scalefactors_json.json'))
spot.size = 16 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data
#[1] 16.0005
#This means that the minimum distance between two neighboring spots is 16 microns.

#Problem with queryKNN if spatial.locs is a data.frame:
spatial.locs <- as.matrix(spatial.locs)


# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

# Set the human ligand-receptor interaction database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Preprocessing
cellchat <- subsetData(cellchat)  # subset to signaling genes
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network
# Compute communication probabilities with spatial constraint
# when using scale.distance=0.01 got an error whith this message: Please increase the value of `scale.distance` and use a value that is slighly smaller than 4.4


cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 50, scale.distance = 4.3, 
                              contact.dependent = TRUE, contact.range = 16)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.table(df.net, file= "Patient_38_CellChat_Interactions_final_annotation.txt", sep='\t')

# Compute communication at pathway level and aggregate
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# Visualization

# Heatmap of communication
pdf("Patient_38_CellChat_Heatmap_final_annotation.pdf")
netVisual_heatmap(cellchat)
dev.off()

# Circle plot showing overall communication between compartments
pdf("Patient_38_CellChat_Circle_final_annotation.pdf")
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, remove.isolate = TRUE)
dev.off()

idents_Post <- grep("_Post$", levels(cellchat@meta$labels), value = TRUE)
idents_Pre <- grep("_Pre$", levels(cellchat@meta$labels), value = TRUE)

pdf("Patient_38_CellChat_Circle_Pre_final_annotation.pdf")
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, idents.use = idents_Pre, remove.isolate = TRUE)
dev.off()

pdf("Patient_38_CellChat_Circle_Post_final_annotation.pdf")
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, idents.use = idents_Post, remove.isolate = TRUE)
dev.off()


# Save CellChat object
saveRDS(cellchat, file = "Patient_38_CellChat_final_annotation.rds")



#Sorting cellchat levels to improve plot visualization

idents <- levels(cellchat@idents)

pre  <- sort(idents[grepl("_Pre$", idents)])
post <- sort(idents[grepl("_Post$", idents)])

new_idents <- c(pre, post)

cellchat@idents <- factor(cellchat@idents, levels = new_idents)

levels(cellchat@idents)


# Final filtered plots

idents_to_filter <- c("B_cells_Pre", "Endothelial_cells_Pre", "Myeloid_Pre",  
                      "Stromal_Pre",  "T_cells_Pre","tumor_Pre",
                      "B_cells_Post", "Endothelial_cells_Post", "Myeloid_Post",
                      "Stromal_Post", "T_cells_Post", "tumor_Post")

pre_idents  <- sort(idents_to_filter[grepl("_Pre$", idents_to_filter)])
post_idents <- sort(idents_to_filter[grepl("_Post$", idents_to_filter)])

# Heatmap

png("Patient_38_cellchat_heatmap.png", width = 4500, height = 4500, res = 500)
netVisual_heatmap(
  cellchat,
  measure = "weight",
  row.show = idents_to_filter,
  col.show = idents_to_filter,
  remove.isolate = TRUE,
  font.size = 18,
  font.size.title = 22,
  ) 
dev.off()

png("Patient_38_CellChat_Circle_Pre_filtered.png", width = 4500, height = 4500, res = 800)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = pre_idents, sources.use = pre_idents, remove.isolate = TRUE, vertex.label.cex = 1.5,
                 edge.label.cex = 1.5, arrow.size = 0.5)
dev.off()

png("Patient_38_CellChat_Circle_Post_filtered.png", width = 4500, height = 4500, res = 800)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge = FALSE, 
                 targets.use = post_idents, sources.use = post_idents, remove.isolate = TRUE, vertex.label.cex = 1.5,
                 edge.label.cex = 1.5, arrow.size = 0.5)
dev.off()

