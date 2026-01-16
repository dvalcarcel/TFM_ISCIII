library(semla)
library(tibble)
library(ggplot2)
library(patchwork)
library(scico)
library(tidyr)
library(dplyr)
library(Seurat)
library(stringr)


setwd("PATH")

load("Patient_38_RCTD_snRNAseq.rda")

# To split object by Denosumab_cluster annotation
objs_list <- SplitObject(Patient_38_RCTD, split.by = "Denosumab_cluster")
Pre1_38  <- objs_list[["Pre1_38"]]
Pre2_38 <- objs_list[["Pre2_38"]]
Post_38 <- objs_list[["Post_38"]]

# Checking the splitted objects
SpatialDimPlot(Pre1_38, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
SpatialDimPlot(Pre2_38, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
SpatialDimPlot(Post_38, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)


# Saving
saveRDS(Pre1_38, "Patient38_Pre1.rds")
saveRDS(Pre2_38, "Patient38_Pre2.rds")
saveRDS(Post_38, "Patient38_Post.rds")


# To plot the Denosumab_cluster annotation to understand how will radial distances be performed

Patient_38_RCTD$Denosumab_cluster[
  !Patient_38_RCTD$Denosumab_cluster %in% c("Pre1_38", "Pre2_38", "Post_38")
] <- NA


png("Patient_38_SpatialDimPlot_Denosumab_cluster.png", width = 4000, height = 1200, res = 450)

SpatialDimPlot(
  Patient_38_RCTD,
  group.by = "Denosumab_cluster",
  image.alpha = 0,
  pt.size.factor = 3.2
) +
  ggtitle("A) Patient 38 – Denosumab cluster annotation") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),  
      ncol = 2
    )
  )
dev.off()

######### Patient 38 pre 1 #############################

Pre1_38 <- Pre1_38[, !is.na(Pre1_38$final_annotation)]

Idents(Pre1_38) <- Pre1_38$final_annotation

Pre1_38_semla <- UpdateSeuratForSemla(Pre1_38)

Pre1_38_semla <- LoadImages(Pre1_38_semla)


pdf("Patient_38_pre1_MapLabels.pdf", height = 8, width = 15)
MapLabels(Pre1_38_semla, column_name = "final_annotation", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

Pre1_38_semla$tumor <- ifelse(Pre1_38_semla$final_annotation %in% "tumor", "tumor", NA)
pdf("Patient_38_pre1_MapLabels_tumor.pdf", height = 8, width = 15)
MapLabels(Pre1_38_semla, column_name = "tumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

# Radial distances from tumor spots

Pre1_38_semla <- RadialDistance(Pre1_38_semla, column_name = "final_annotation", selected_groups = "tumor")

pdf("Patient_38_pre1_MapFeatures_tumor.pdf", height = 8, width = 15)
MapFeatures(Pre1_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()

#Pixel coordinates 

Pre1_38_semla$r_dist_tumor <- (100/273)*Pre1_38_semla$r_dist_tumor

Pre1_38_semla <- RadialDistance(Pre1_38_semla, column_name = "final_annotation", 
                     selected_groups = "tumor", convert_to_microns = TRUE)


Pre1_38_semla$r_dist_tumor_sqrt <- sign(Pre1_38_semla$r_dist_tumor)*sqrt(abs(Pre1_38_semla$r_dist_tumor))
pdf("Patient_38_pre1_MapFeatures_tumor_2.pdf", height = 8, width = 15)
MapFeatures(Pre1_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()


saveRDS(Pre1_38_semla, "Patient38_Pre1_semla.rds")



######### Patient 38 pre 2 #############################

Pre2_38 <- Pre2_38[, !is.na(Pre2_38$final_annotation)]

Idents(Pre2_38) <- Pre2_38$final_annotation

Pre2_38_semla <- UpdateSeuratForSemla(Pre2_38)

Pre2_38_semla <- LoadImages(Pre2_38_semla)


pdf("Patient_38_pre2_MapLabels.pdf", height = 8, width = 15)
MapLabels(Pre2_38_semla, column_name = "final_annotation", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

Pre2_38_semla$tumor <- ifelse(Pre2_38_semla$final_annotation %in% "tumor", "tumor", NA)
pdf("Patient_38_pre2_MapLabels_tumor.pdf", height = 8, width = 15)
MapLabels(Pre2_38_semla, column_name = "tumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()


# Radial distances from tumor spots

Pre2_38_semla <- RadialDistance(Pre2_38_semla, column_name = "final_annotation", selected_groups = "tumor")

pdf("Patient_38_pre2_MapFeatures_tumor.pdf", height = 8, width = 15)
MapFeatures(Pre2_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()

#Pixel coordinates 

Pre2_38_semla$r_dist_tumor <- (100/273)*Pre2_38_semla$r_dist_tumor

Pre2_38_semla <- RadialDistance(Pre2_38_semla, column_name = "final_annotation", 
                     selected_groups = "tumor", convert_to_microns = TRUE)


Pre2_38_semla$r_dist_tumor_sqrt <- sign(Pre2_38_semla$r_dist_tumor)*sqrt(abs(Pre2_38_semla$r_dist_tumor))
pdf("Patient_38_pre2_MapFeatures_tumor_2.pdf", height = 8, width = 15)
MapFeatures(Pre2_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()

saveRDS(Pre2_38_semla, "Patient38_Pre2_semla.rds")


######### Patient 38 post #############################

Post_38 <- Post_38[, !is.na(Post_38$final_annotation)]

Idents(Post_38) <- Post_38$final_annotation

Post_38_semla <- UpdateSeuratForSemla(Post_38)

Post_38_semla <- LoadImages(Post_38_semla)


pdf("Patient_38_post_MapLabels.pdf", height = 8, width = 15)
MapLabels(Post_38_semla, column_name = "final_annotation", 
          image_use = "raw", pt_alpha = 0.6, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()

Post_38_semla$tumor <- ifelse(Post_38_semla$final_annotation %in% "tumor", "tumor", NA)
pdf("Patient_38_post_MapLabels_tumor.pdf", height = 8, width = 15)
MapLabels(Post_38_semla, column_name = "tumor", override_plot_dims = TRUE, 
          image_use = "raw", drop_na = TRUE, pt_size = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  guides(fill = guide_legend(override.aes = list(size = 3), ncol = 2))
dev.off()


# Radial distances from tumor spots

Post_38_semla <- RadialDistance(Post_38_semla, column_name = "final_annotation", selected_groups = "tumor")

pdf("Patient_38_post_MapFeatures_tumor.pdf", height = 8, width = 15)
MapFeatures(Post_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()

#Pixel coordinates 

Post_38_semla$r_dist_tumor <- (100/273)*Post_38_semla$r_dist_tumor

Post_38_semla <- RadialDistance(Post_38_semla, column_name = "final_annotation", 
                                 selected_groups = "tumor", convert_to_microns = TRUE)


Post_38_semla$r_dist_tumor_sqrt <- sign(Post_38_semla$r_dist_tumor)*sqrt(abs(Post_38_semla$r_dist_tumor))
pdf("Patient_38_post_MapFeatures_tumor_2.pdf", height = 8, width = 15)
MapFeatures(Post_38_semla, features = "r_dist_tumor", center_zero = TRUE, pt_size = 2, 
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev(),
            override_plot_dims = TRUE)
dev.off()


#### Radial distances Tumor - Myeloid, T Cells, B cells and Endothelial cells
Pre1_38_semla <- readRDS(file = "Patient38_Pre1_semla.rds")
Pre2_38_semla <- readRDS(file = "Patient38_Pre2_semla.rds")
Post_38_semla <- readRDS(file = "Patient38_Post_semla.rds")


df_merged_38 <- bind_rows(
  Pre1_38_semla@meta.data |> 
    as_tibble(rownames = "barcode") |>  
    mutate(sample = "pre1"),
  
  Pre2_38_semla@meta.data |> 
    as_tibble(rownames = "barcode") |>  
    mutate(sample = "pre2"),
  
  Post_38_semla@meta.data |> 
    as_tibble(rownames = "barcode") |>  
    mutate(sample = "post")
)

saveRDS(df_merged_38, "df_merged_38_semla.rds")


sample_colors <- c(
  "pre1"           = "#2ca02c",
  "pre2"      = "#1f78b4",
  "post"       = "#d62728")


# FIltering, pivoting and plotting

sel_pathway <- c("T cells") 

p1 <- df_merged_38 |>  
  select(barcode, sample, r_dist_tumor, all_of(sel_pathway)) |>  
  filter(r_dist_tumor < 1e3) |>  
  drop_na(all_of(sel_pathway)) |>  
  pivot_longer(cols = all_of(sel_pathway), names_to = "cell_type", values_to = "score") |>  
  ggplot(aes(x = r_dist_tumor, y = score, color = sample)) +   
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +  
  geom_vline(aes(xintercept = 0), linetype = "dashed") +  
  theme_minimal() + 
  scale_color_manual(values = sample_colors) +
  labs(x = "Radial Distance Tumor (µm)", 
       y = "Pathway Score", 
       title = "B) Patient 38 - Tumor Radial Distances - T cells" ,
       color = "Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(limits = c(-1000, 1000)) +
  scale_y_continuous(limits = c(-0.02, 0.4))

sel_pathway <- c("Myeloid") 

p2 <- df_merged_38 |>  
  select(barcode, sample, r_dist_tumor, all_of(sel_pathway)) |>  
  filter(r_dist_tumor < 1e3) |>  
  drop_na(all_of(sel_pathway)) |>  
  pivot_longer(cols = all_of(sel_pathway), names_to = "cell_type", values_to = "score") |>  
  ggplot(aes(x = r_dist_tumor, y = score, color = sample)) +   
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +  
  geom_vline(aes(xintercept = 0), linetype = "dashed") +  
  theme_minimal() + 
  scale_color_manual(values = sample_colors) +
  labs(x = "Radial Distance Tumor (µm)", 
       y = "Pathway Score", 
       title = "C) Patient 38 - Tumor Radial Distances - Myeloid" ,
       color = "Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(limits = c(-1000, 1000)) +
  scale_y_continuous(limits = c(-0.02, 0.4))



sel_pathway <- c("B cells") 

p3 <- df_merged_38 |>  
  select(barcode, sample, r_dist_tumor, all_of(sel_pathway)) |>  
  filter(r_dist_tumor < 1e3) |>  
  drop_na(all_of(sel_pathway)) |>  
  pivot_longer(cols = all_of(sel_pathway), names_to = "cell_type", values_to = "score") |>  
  ggplot(aes(x = r_dist_tumor, y = score, color = sample)) +   
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +  
  geom_vline(aes(xintercept = 0), linetype = "dashed") +  
  theme_minimal() + 
  scale_color_manual(values = sample_colors) +
  labs(x = "Radial Distance Tumor (µm)", 
       y = "Pathway Score", 
       title = "D) Patient 38 - Tumor Radial Distances - B cells" ,
       color = "Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(limits = c(-1000, 1000)) +
  scale_y_continuous(limits = c(-0.02, 0.4))


sel_pathway <- c("Endothelial cells") 

p4 <- df_merged_38 |>  
  select(barcode, sample, r_dist_tumor, all_of(sel_pathway)) |>  
  filter(r_dist_tumor < 1e3) |>  
  drop_na(all_of(sel_pathway)) |>  
  pivot_longer(cols = all_of(sel_pathway), names_to = "cell_type", values_to = "score") |>  
  ggplot(aes(x = r_dist_tumor, y = score, color = sample)) +   
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +  
  geom_vline(aes(xintercept = 0), linetype = "dashed") +  
  theme_minimal() + 
  scale_color_manual(values = sample_colors) +
  labs(x = "Radial Distance Tumor (µm)", 
       y = "Pathway Score", 
       title = "E) Patient 38 - Tumor Radial Distances - Endothelial cells" ,
       color = "Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(limits = c(-1000, 1000)) +
  scale_y_continuous(limits = c(-0.02, 0.4))




png("Patient_38_radial_distances_tumor.png", width = 6000, height = 3600, res = 450)

p1+p2+p3+p4

dev.off()


############## Interferon GAMMA pre and post #######

sel_pathway <- c("ONLY_EXP_UP_IFN_GAMMA") 

png("Patient_38_merged_tumor_IFN_GAMMA.png", width = 6000, height = 2400, res = 450)

df_merged_38 |>  
  select(barcode, sample, r_dist_tumor, all_of(sel_pathway)) |>  
  filter(r_dist_tumor < 1e3) |>  
  drop_na(all_of(sel_pathway)) |>  
  pivot_longer(cols = all_of(sel_pathway), names_to = "cell_type", values_to = "score") |>  
  ggplot(aes(x = r_dist_tumor, y = score, color = sample)) +   
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +  
  geom_vline(aes(xintercept = 0), linetype = "dashed") +  
  theme_minimal() + 
  scale_color_manual(values = sample_colors) +
  labs(x = "Radial Distance Tumor (µm)", 
       y = "Pathway Score", 
       title = "F) Patient 38 - Tumor Radial Distances - HALLMARK_INTERFERON_GAMMA_RESPONSE (Bulk-RNAseq)" ,
       color = "Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  scale_x_continuous(limits = c(-1000, 1000)) +
  scale_y_continuous(limits = c(-0.02, 0.4))
dev.off()



