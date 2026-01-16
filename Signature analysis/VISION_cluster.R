
setwd("CLUSTER PATH")

library("Seurat")
library("VISION")
library("ggplot2")

# This script was run on the CNIO's cluster:


load("Patient_38_RCTD_snRNAseq.rda")

Patient_38_RCTD_projection <- Patient_38_RCTD@reductions$umap@cell.embeddings

# We need to previously create a gmt file with the signatures

Patient_38_RCTD_vis <- Vision(Patient_38_RCTD@assays$Spatial$data,signatures=c("GMT file"), meta = Patient_38_RCTD@meta.data, pool=FALSE)
Patient_38_RCTD_vis <- addProjection(Patient_38_RCTD_vis, "UMAP", Patient_38_RCTD_projection)
Patient_38_RCTD_vis <- analyze(Patient_38_RCTD_vis)

save("Patient_38_RCTD_vis", file = "Patient_38_RCTD_vis.rda")
