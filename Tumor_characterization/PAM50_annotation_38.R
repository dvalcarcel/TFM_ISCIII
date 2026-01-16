if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefu")

library("genefu")
data(pam50)
data(pam50.robust)
library("Seurat")
library("org.Hs.eg.db")
library("dplyr")
library("ggplot2")
library("gghalves")
library("tidyverse")
library("patchwork")


setwd("PATH")

load("Patient_38_RCTD_snRNAseq.rda")

# We need to select only tumor spots

tumor_spots <- WhichCells(Patient_38_RCTD, ident = "tumor")

mat <- GetAssayData(Patient_38_RCTD, assay = "Spatial", slot = "data")[, tumor_spots]

mat_subtype <- t(as.matrix(mat))  # spots × genes

# We need a matrix with the spots as rownames and genes as colnames
dim(mat_subtype)
head(rownames(mat_subtype))   # spots
head(colnames(mat_subtype))   # genes

genes <- colnames(mat_subtype)

entrez <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

annot <- data.frame(
  probe = genes,
  EntrezGene.ID = entrez,
  row.names = genes
)

head(annot)
table(is.na(annot$EntrezGene.ID))

res <- molecular.subtyping(
  sbt.model = "pam50",
  data = mat_subtype,
  annot = annot,
  do.mapping = TRUE,
  verbose = TRUE
)



# Now we need to annotate the probability of each PAM50 subtype in order to 
# have an annotation for all the subtypes

probs <- res$subtype.proba


for (subtype in colnames(probs)) {
  
  Patient_38_RCTD[[paste0("PAM50_prob_", subtype)]] <- NA
  
  Patient_38_RCTD@meta.data[tumor_spots, paste0("PAM50_prob_", subtype)] <- 
    probs[tumor_spots, subtype]
}


# We also need to add a new annotation with the PAM50 dominant subtype per spot

pam50_dominant <- as.data.frame(res$subtype)
colnames(pam50_dominant) <- "PAM50_dominant_subtype"
Patient_38_RCTD <- AddMetaData(Patient_38_RCTD, pam50_dominant)

# We need to check if the PAM50 annotation within the VisiumHD object is equal to
# the genefu annotation

print(table(Patient_38_RCTD$PAM50_dominant_subtype, useNA = "ifany"))
print(table(res$subtype))

save(Patient_38_RCTD, file = "Patient_38_RCTD_snRNAseq.rda")


# Define output directory for signature plots
output_dir <- "signature_plots_38/"

# To define signature names 
signature_names <- c("PAM50_prob_Basal", "PAM50_prob_Her2", "PAM50_prob_LumA",
                     "PAM50_prob_LumB", "PAM50_prob_Normal")

# Generate and save SpatialFeaturePlots for each signature
for (signature in signature_names) {
  p <- SpatialFeaturePlot(Patient_38_RCTD, features = signature, image.alpha=0,pt.size.factor=3.2) + ggtitle(signature)
  
  # Save the plot as a PDF
  pdf_filename <- paste0(output_dir, signature, ".pdf")
  ggsave(pdf_filename, plot = p, width = 10, height = 7)
  
  print(paste("Saved:", pdf_filename))
}

#SpatialDimPlot to spatially visualize the PAM50 dominant subtype

p0 <- SpatialDimPlot(Patient_38_RCTD, group.by ="PAM50_dominant_subtype", image.alpha=0, pt.size.factor=3.2) +
  ggtitle("A) Patient 38 – SpatialDimPlot - PAM50 dominant subtype") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text  = element_text(size = 22)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 16),  
      ncol = 2
    ))


# We are going to visualize each PAM50 subtype score through a VlnPlot and their
# comparison with the treatment

obj_tumor <- subset(Patient_38_RCTD, idents = "tumor")

df <- obj_tumor@meta.data

df_long <- df %>%
  pivot_longer(
    cols = c(
      PAM50_prob_Basal,
      PAM50_prob_Her2,
      PAM50_prob_LumA,
      PAM50_prob_LumB,
      PAM50_prob_Normal
    ),
    names_to = "Score",
    values_to = "Value"
  )

df_long$Denosumab <- factor(
  df_long$Denosumab,
  levels = c("Pre", "Post")
)

p1 <- ggplot(df_long, aes(x = Patient, y = Value, fill = Denosumab)) +
  geom_half_violin(
    data = subset(df_long, Denosumab == "Pre"),
    side = "l",
    alpha = 0.6
  ) +
  geom_half_violin(
    data = subset(df_long, Denosumab == "Post"),
    side = "r",
    alpha = 0.6
  ) +
  facet_wrap(~Score, scales = "free_y") +
  scale_fill_manual(values = c("Pre" = "#1f77b4", "Post" = "#ff7f0e")) +
  theme_minimal() +
  ylab("Score") +
  xlab("Patient") +
  ggtitle("B) Patient 38 - PAM50 scores: Pre vs Post Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20),
    axis.text.x  = element_text(size = 20),  
    axis.text.y  = element_text(size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold")
  )



# Finally, we are going to visualize each PAM50 subtype score through a barplot and their
# comparison with the treatment, normalizing by the proportion of spots

df$Denosumab <- factor(
  df$Denosumab,
  levels = c("Pre", "Post")
)

Patient_38_RCTD$final_annotation_Denosumab <-paste(Patient_38_RCTD$final_annotation, Patient_38_RCTD$Denosumab, sep = "_")


df_tumor <- Patient_38_RCTD@meta.data %>%
  filter(grepl("^tumor_P", final_annotation_Denosumab))

df_tumor <- df_tumor %>%
  mutate(
    Denosumab = case_when(
      grepl("_Pre$", final_annotation_Denosumab)  ~ "Pre",
      grepl("_Post$", final_annotation_Denosumab) ~ "Post",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Denosumab))


df_prop <- df_tumor %>%
  count(Denosumab, PAM50_dominant_subtype) %>%
  group_by(Denosumab) %>%
  mutate(
    proportion = n / sum(n)
  )

df_prop$Denosumab <- factor(
  df_prop$Denosumab,
  levels = c("Pre", "Post")
)


p2 <- ggplot(df_prop,
             aes(x = PAM50_dominant_subtype, y = proportion, fill = Denosumab)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Pre" = "#1f77b4", "Post" = "#ff7f0e")) +
  theme_minimal() +
  ylab("Proportion (normalized by number of tumor spots)") +
  xlab("PAM50 dominant subtype") +
  ggtitle("C) Patient 38 - Dominant PAM50 subtypes normalized by number of tumor spots\nPre (n=4445) vs Post (n=13340)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text  = element_text(size = 22),
    axis.text.x  = element_text(size = 22),  
    axis.text.y  = element_text(size = 22),  
    axis.title.x = element_text(size = 22),  
    axis.title.y = element_text(size = 22)
  )




tab <- table(df_tumor$Denosumab, df_tumor$PAM50_dominant_subtype)
chisq.test(tab)


p_all <- (p0) / (p1 | p2) 


pdf("Patient_38_PAM50_subtypes.pdf", width = 28, height = 32)
p_all
dev.off()


png("Patient_38_PAM50_subtypes.png", width = 28, height = 32, units = "in", res = 600)
p_all
dev.off()
