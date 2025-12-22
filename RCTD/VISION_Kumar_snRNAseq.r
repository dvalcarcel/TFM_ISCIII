
setwd("PATH/")

library("Seurat")
library("VISION")
library("ggplot2")

load("snRNAseq_merged_integrated.rda")


snRNAseq_projection <- snRNAseq_merged_integrated@reductions$umap@cell.embeddings
snRNAseq_vis <- Vision(snRNAseq_merged_integrated@assays$RNA$data,signatures=c("./sn_top_markers_Kumar_et_al.gmt"), meta = snRNAseq_merged_integrated@meta.data, pool=FALSE)
snRNAseq_vis <- addProjection(snRNAseq_vis, "UMAP", snRNAseq_projection)
snRNAseq_vis <- analyze(snRNAseq_vis)

save("snRNAseq_vis", file = "snRNAseq_vis.rda")

signature_scores <- as.data.frame(snRNAseq_vis@SigScores)

snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Basal, col.name = "Basal")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$LumHR, col.name = "LumHR")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$LumSec, col.name = "LumSec")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Fibroblasts, col.name = "Fibroblasts")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Lymphatic, col.name = "Lymphatic")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Mast, col.name = "Mast")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Vascular, col.name = "Vascular")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Perivascular, col.name = "Perivascular")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Adipocytes, col.name = "Adipocytes")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$T_cells, col.name = "T_cells")
snRNAseq_merged_integrated <- AddMetaData(snRNAseq_merged_integrated, metadata = signature_scores$Myeliod, col.name = "Myeliod")

basal <- FeaturePlot(snRNAseq_merged_integrated, features = "Basal", reduction = "umap")
pdf("snRNAseq_Basal.pdf")
basal + scale_fill_gradientn(colors = c("blue", "white", "red"))  
dev.off() 
lumHR <- FeaturePlot(snRNAseq_merged_integrated, features = "LumHR", reduction = "umap")
pdf("snRNAseq_LumHR.pdf")
lumHR + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
lumSec <- FeaturePlot(snRNAseq_merged_integrated, features = "LumSec", reduction = "umap")
pdf("snRNAseq_LumSec.pdf")
lumSec + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
fibroblasts <- FeaturePlot(snRNAseq_merged_integrated, features = "Fibroblasts", reduction = "umap")
pdf("snRNAseq_Fibroblasts.pdf")
fibroblasts + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
lymphatic <- FeaturePlot(snRNAseq_merged_integrated, features = "Lymphatic", reduction = "umap")
pdf("snRNAseq_Lymphatic.pdf")
lymphatic + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
mast <- FeaturePlot(snRNAseq_merged_integrated, features = "Mast", reduction = "umap")
pdf("snRNAseq_Mast.pdf")
mast + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
vascular <- FeaturePlot(snRNAseq_merged_integrated, features = "Vascular", reduction = "umap")
pdf("snRNAseq_Vascular.pdf")
vascular + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
perivascular <- FeaturePlot(snRNAseq_merged_integrated, features = "Perivascular", reduction = "umap")
pdf("snRNAseq_Perivascular.pdf")
perivascular + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
adipocytes <- FeaturePlot(snRNAseq_merged_integrated, features = "Adipocytes", reduction = "umap")
pdf("snRNAseq_Adipocytes.pdf")
adipocytes + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
t_cells <- FeaturePlot(snRNAseq_merged_integrated, features = "T_cells", reduction = "umap")
pdf("snRNAseq_T_cells.pdf")
t_cells + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()
myeliod <- FeaturePlot(snRNAseq_merged_integrated, features = "Myeliod", reduction = "umap")
pdf("snRNAseq_Myeliod.pdf")
myeliod + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

save("snRNAseq_merged_integrated", file = "snRNAseq_merged_integrated.rda")






