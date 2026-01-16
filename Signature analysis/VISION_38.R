library("Seurat")
library("VISION")
library("ggplot2")

setwd("PATH")

load("Patient_38_RCTD_snRNAseq.rda")

#We need to load the VISION object created from the "VISION_cluster" script:
load("Patient_38_RCTD_vis.rda")


signature_scores <- as.data.frame(Patient_38_RCTD_vis@SigScores)

all(rownames(signature_scores) %in% colnames(Patient_38_RCTD))

Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, metadata = signature_scores)


# Define output directory for PDFs
output_dir <- "signature_plots_38/"
dir.create(output_dir, showWarnings = FALSE)

# Extract signature names from metadata
signature_names <- colnames(signature_scores)

# Generate and save SpatialFeaturePlots for each signature
for (signature in signature_names) {
  p <- SpatialFeaturePlot(Patient_38_RCTD, features = signature,image.alpha=0,pt.size.factor=3.2) + ggtitle(signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, ".pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}

save(Patient_38_RCTD, file = "Patient_38_RCTD_snRNAseq.rda")


# Generate and save VlnPlots for each signature
for (signature in signature_names) {
  p <- VlnPlot(Patient_38_RCTD, features = signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, "_VlnPlot.pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}


