library(stats)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(patchwork)


#To make barplot function for GSEA's results, applying a FDR threshold of 5% due to permutation by gene set.
barplot_gsea <- function(df, nes_col = "NES", fdr_col = "FDR q-val", 
                         name_col = "NAME", fdr_threshold = 0.05,
                         title = "GSEA barplot", subtitle = NULL) {
  
  # Rename columns, FDR filtering and sortering by FDR
  df_clean <- df %>%
    rename(NES = all_of(nes_col),
           FDR = all_of(fdr_col),
           NAME = all_of(name_col)) %>%
    filter(FDR < fdr_threshold) %>%
    arrange(FDR)
  
  # Plot
  p <- ggplot(df_clean, aes(x = NES, y = reorder(NAME, NES), fill = FDR)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue") + 
    labs(fill = "FDR q-val",
         title = title,
         x = "Normalized Enrichment Score (NES)",
         y = "Pathway") +
    scale_x_continuous(limits = c(-3, 3)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.text  = element_text(size = 20),
          axis.text.x  = element_text(size = 20),  
          axis.text.y  = element_text(size = 20),  
          axis.title.x = element_text(size = 20),  
          axis.title.y = element_text(size = 20)) +
    annotate("text", x = 0, y = -0.5, 
             label = "Pre (Negative NES values) | Post (Positive NES values)",
             vjust = -1.5, size = 7, color = "gray30", fontface = "italic")
  
  return(p)
}


# Set path with the outputs of GSEA
setwd("PATH")

#Tumor cells hallmark
tumor_cells_hall_post <- read_tsv(file = "PB_Tumor_cells.Gsea.1761724193240/gsea_report_for_tumor_post_1761724193240.tsv")
tumor_cells_hall_pre <- read_tsv(file = "PB_Tumor_cells.Gsea.1761724193240/gsea_report_for_tumor_pre_1761724193240.tsv")
tumor_cells_hall_post$Condition <- "Post"
tumor_cells_hall_pre$Condition  <- "Pre"

tumor_cells_hall_merged <- bind_rows(tumor_cells_hall_pre, tumor_cells_hall_post)


tumor_cells_hall_barplot <- barplot_gsea(
  tumor_cells_hall_merged,
  title = "B) GSEA barplot - Tumor cells Pre vs Post (Hallmark)"
)
png("GSEA_Tumor_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
tumor_cells_hall_barplot
dev.off()


#T cells hallmark
T_cells_hall_post <- read_tsv(file = "PB_T_cells.Gsea.1761723797030/gsea_report_for_T_post_1761723797030.tsv")
T_cells_hall_pre <- read_tsv(file = "/Users/dvalcarcel/Documents/TFM/GSEA_PB/oct28/PB_T_cells.Gsea.1761723797030/gsea_report_for_T_pre_1761723797030.tsv")
T_cells_hall_post$Condition <- "Post"
T_cells_hall_pre$Condition  <- "Pre"

T_cells_hall_merged <- bind_rows(T_cells_hall_pre, T_cells_hall_post)


T_cells_hall_barplot <- barplot_gsea(
  T_cells_hall_merged,
  title = "B) GSEA barplot - T cells Pre vs Post (Hallmark)")
png("GSEA_T_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
T_cells_hall_barplot
dev.off()


#Myeloid cells hallmark
Myeloid_cells_hall_post <- read_tsv(file = "PB_Myeloid_cells.Gsea.1761666543638/gsea_report_for_Myeloid_post_1761666543638.tsv")
Myeloid_cells_hall_pre <- read_tsv(file = "PB_Myeloid_cells.Gsea.1761666543638/gsea_report_for_Myeloid_pre_1761666543638.tsv")
Myeloid_cells_hall_post$Condition <- "Post"
Myeloid_cells_hall_pre$Condition  <- "Pre"
Myeloid_cells_hall_merged <- bind_rows(Myeloid_cells_hall_pre, Myeloid_cells_hall_post)


Myeloid_cells_hall_barplot <- barplot_gsea(
  Myeloid_cells_hall_merged,
  title = "D) GSEA barplot - Myeloid cells Pre vs Post (Hallmark)")

png("GSEA_Myeloid_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
Myeloid_cells_hall_barplot
dev.off()

png("DESeq2_GSEA_Myeloid_Post_vs_Pre.png", width = 20000, height = 8000, res = 750)
(p3 | Myeloid_cells_hall_barplot)
dev.off()

#Endothelial cells hallmark
Endothelial_cells_hall_post <- read_tsv(file = "PB_Endothelial_cells.Gsea.1761665389252/gsea_report_for_Endothelial_post_1761665389252.tsv")
Endothelial_cells_hall_pre <- read_tsv(file = "PB_Endothelial_cells.Gsea.1761665389252/gsea_report_for_Endothelial_pre_1761665389252.tsv")
Endothelial_cells_hall_post$Condition <- "Post"
Endothelial_cells_hall_pre$Condition  <- "Pre"
Endothelial_cells_hall_merged <- bind_rows(Endothelial_cells_hall_pre, Endothelial_cells_hall_post)


Endothelial_cells_hall_barplot <- barplot_gsea(
  Endothelial_cells_hall_merged,
  title = "F) GSEA barplot - Endothelial cells Pre vs Post (Hallmark)")

png("GSEA_Endothelial_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
Endothelial_cells_hall_barplot
dev.off()

png("DESeq2_GSEA_Endothelial_Post_vs_Pre.png", width = 20000, height = 8000, res = 750)
(p4 | Endothelial_cells_hall_barplot)
dev.off()


#B cells hallmark
B_cells_hall_post <- read_tsv(file = "PB_B_cells.Gsea.1761665074107/gsea_report_for_B_post_1761665074107.tsv")
B_cells_hall_pre <- read_tsv(file = "PB_B_cells.Gsea.1761665074107/gsea_report_for_B_pre_1761665074107.tsv")
B_cells_hall_post$Condition <- "Post"
B_cells_hall_pre$Condition  <- "Pre"
B_cells_hall_merged <- bind_rows(B_cells_hall_pre, B_cells_hall_post)


B_cells_hall_barplot <- barplot_gsea(
  B_cells_hall_merged,
  title = "D) GSEA barplot - B cells Pre vs Post (Hallmark)")

png("GSEA_B_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
B_cells_hall_barplot
dev.off()


# Supplementary figure with DESeq2 and GSEA plots for T and B cells
p_all <- (p1 | T_cells_hall_barplot) / (p2 | B_cells_hall_barplot)
pdf("DESeq2_GSEA_T_B_cells.pdf", width = 28, height = 32)
p_all
dev.off()



#Stromal cells hallmark
Stromal_cells_hall_post <- read_tsv(file = "PB_Stromal_cells.Gsea.1761723579756/gsea_report_for_Stromal_post_1761723579756.tsv")
Stromal_cells_hall_pre <- read_tsv(file = "PB_Stromal_cells.Gsea.1761723579756/gsea_report_for_Stromal_pre_1761723579756.tsv")
Stromal_cells_hall_post$Condition <- "Post"
Stromal_cells_hall_pre$Condition  <- "Pre"
Stromal_cells_hall_merged <- bind_rows(Stromal_cells_hall_pre, Stromal_cells_hall_post)


Stromal_cells_hall_barplot <- barplot_gsea(
  Stromal_cells_hall_merged,
  title = "H) GSEA barplot - Stromal cells Pre vs Post (Hallmark)")

png("GSEA_Stromal_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)
Stromal_cells_hall_barplot
dev.off()

png("DESeq2_GSEA_Stromal_Post_vs_Pre.png", width = 20000, height = 8000, res = 750)
(p5 | Stromal_cells_hall_barplot)
dev.off()


