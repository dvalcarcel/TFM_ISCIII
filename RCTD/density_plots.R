library(ggplot2)
library(ggridges)
library(tidyr)
library(dplyr)
library(viridis)


setwd("PATH")

############ Patient 38 ###################################

load("Patient_38_RCTD_snRNAseq.rda")

epithelial.df <- data.frame(
  Epithelial_cells = Patient_38_RCTD$`Epithelial cells`,
  Malignant = Patient_38_RCTD$malignant,
  Basal = Patient_38_RCTD$Basal,
  LumSec = Patient_38_RCTD$LumSec,
  LumHR = Patient_38_RCTD$LumHR,
  Pathologist_annotation_refined = Patient_38_RCTD$Pathologist_annotation_refined
)

Idents(Patient_38_RCTD) <- Patient_38_RCTD$Pathologist_annotation_refined

epithelial.df %>%
  group_by(Pathologist_annotation_refined) %>%
  dplyr::filter(Pathologist_annotation_refined == "tumor") %>%
  ungroup() -> tumor_patient_38

epithelial.df %>%
  group_by(Pathologist_annotation_refined) %>%
  dplyr::filter(Pathologist_annotation_refined == "normal_breast_tissue") %>%
  ungroup() -> normal_patient_38


# Convert to long format
tumor_epithelial_long <- tumor_patient_38 %>%
  pivot_longer(cols = c(Epithelial_cells, Malignant, Basal, LumSec, LumHR), 
               names_to = "Cell_Type", 
               values_to = "Weight")

normal_epithelial_long <- normal_patient_38 %>%
  pivot_longer(cols = c(Epithelial_cells, Malignant, Basal, LumSec, LumHR), 
               names_to = "Cell_Type", 
               values_to = "Weight")

# Density plots

p1 <- ggplot(tumor_epithelial_long, aes(x = Weight, y = Cell_Type, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Weight", option = "C") +
  labs(
    title = "A) Distribution of Epithelial Cell Types in spots annotated as tumors by the pathologist",
    x = "Weight",
    y = "Cell Type"
  ) +
  theme_minimal()

p2 <- ggplot(normal_epithelial_long, aes(x = Weight, y = Cell_Type, fill = stat(x))) +
   geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
   scale_fill_viridis_c(name = "Weight", option = "C") +
   labs(
     title = "B) Distribution of Epithelial Cell Types in spots annotated as normal breast by the pathologist",
     x = "Weight",
     y = "Cell Type"
   ) +
   theme_minimal()

pdf("Patient_38_Epithelial_density_plot.pdf",  height = 7, width = 20)
p1+p2
dev.off()






