library(dplyr)
library(ggplot2)
library(Hmisc)
library(scales)
library(stats)
library(RColorBrewer)
library(tidyr)
library(Seurat)
library(STDistance)
library(patchwork)

setwd("PATH")


######### Patient 38 #############

load("Patient_38_RCTD_snRNAseq.rda")

# Extract coordinates from the Seurat object
coords_38 <- GetTissueCoordinates(Patient_38_RCTD)

#Create a data frame type “tissue_posi” from the STDistance tutorial
tissue_posi_38 <- data.frame(
  barcode = rownames(coords_38),
  in_tissue = 1,  
  array_row = NA,
  array_col = NA,
  pxl_row_in_fullres = coords_38$y,   
  pxl_col_in_fullres = coords_38$x,   
  Sample = "Patient_38",
  Sampleid = 38,
  Newbarcode = paste0(rownames(coords_38), "_38")
)

# To create metadata from Seurat
meta_38 <- Patient_38_RCTD@meta.data
meta_38$Newbarcode <- paste0(rownames(meta_38), "_38")


# Merging coordinates with metada
posi_38 <- merge(
  x = tissue_posi_38,
  y = meta_38,
  by = "Newbarcode",
  all.y = TRUE
)

posi_38 <- posi_38[!is.na(posi_38$final_annotation), ]


posi_38_Pre <- posi_38 %>% filter(Denosumab %in% "Pre")
posi_38_Post <- posi_38 %>% filter(Denosumab %in% "Post")




########## Tumor cells as center ##################

########## Myeloid cells as targets ###############

# To calculate nearest distances to Myeloid cells

distance_results_38_TumortoMyeloid_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoMyeloid_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "Myeloid",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)


p1 <- plot_distance_boxplot(
  distance_results_38_TumortoMyeloid_Pre,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("A) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to Myeloid cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 80))  


p2 <- plot_distance_boxplot(
  distance_results_38_TumortoMyeloid_Post,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("B) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to Myeloid cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 80))  

png("Patient_38_STDistance_Tumor_Myeloid_boxplot.png", width = 7200, height = 2400, res = 450)
p1+p2
dev.off()


# Radial layout plots

plot_radial_distance(
  distance_results_38_TumortoMyeloid_Pre,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to Myeloid cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


plot_radial_distance(
  distance_results_38_TumortoMyeloid_Post,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to Myeloid cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

# Spatial network plots

visualize_spatial_network(
  posi_38,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Myeloid",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Myeloid" = "#377EB8"),
  alpha = 0.7
) + ggtitle("Patient 38 - Spatial Network from Tumor to Myeloid cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),
      ncol = 2
    ))


visualize_spatial_network(
  posi_38_Pre,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Myeloid",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Myeloid" = "#377EB8"),
  alpha = 0.7
)

visualize_spatial_network(
  posi_38_Post,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Myeloid",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Myeloid" = "#377EB8"),
  alpha = 0.7
)

# Analyzing the differences through t test
t.test(distance_results_38_TumortoMyeloid_Pre$Myeloid, distance_results_38_TumortoMyeloid_Post$Myeloid)



# Spatial Gradient Network IFN GAMMA

p3 <- visualize_spatial_gradient(
  spatial_data = posi_38_Pre,
  sample = "Patient_38",
  gradient_type = "Myeloid",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
) +   scale_y_reverse() +
  ggtitle("C) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in Myeloid cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.text.x  = element_text(size = 16),  
    axis.text.y  = element_text(size = 16),  
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18)
  ) 


p4 <- visualize_spatial_gradient(
  spatial_data = posi_38_Post,
  sample = "Patient_38",
  gradient_type = "Myeloid",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
) +   scale_y_reverse() +
  ggtitle("D) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in Myeloid cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    axis.text.x  = element_text(size = 16),  
    axis.text.y  = element_text(size = 16),  
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18)
  ) 

png("Patient_38_STDistance_Tumor_Myeloid_SpatialGradientNetwork.png", width = 4700, height = 2500, res = 145)
p3+p4
dev.off()



# Correlation plots between distance Tumor to Myeloid and IFN GAMMA expression

merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoMyeloid_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)
merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoMyeloid_Post,
  by = "Newbarcode",
  all.x = FALSE 
)


result_correlation_38_Pre <- calculate_correlations(
  spatial_data = merged_38_IFN_Pre,
  distance_results = distance_results_38_TumortoMyeloid_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Myeloid",
  method = "pearson",
  plot = TRUE,
  plot_title = "E) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) \n Expression and distance from Tumor to Myeloid cells before Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Pre$estimate))
print(paste("P-value:", result_correlation_38_Pre$p_value))

p5 <- result_correlation_38_Pre$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 

result_correlation_38_Post <- calculate_correlations(
  spatial_data = merged_38_IFN_Post,
  distance_results = distance_results_38_TumortoMyeloid_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Myeloid",
  method = "pearson",
  plot = TRUE,
  plot_title = "F) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq)\n Expression and distance from Tumor to Myeloid cells after Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Post$estimate))
print(paste("P-value:", result_correlation_38_Post$p_value))

p6 <- result_correlation_38_Post$plot + 
  theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 

png("Patient_38_STDistance_Tumor_Myeloid_IFN_correlation.png", width = 4500, height = 2000, res = 225)
p5+p6
dev.off()


p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) 


pdf("Patient_38_STDistance_Tumor_Myeloid.pdf", width = 28, height = 32)
p_all
dev.off()




########## T cells as targets ###############

# To calculate nearest distances to T cells

distance_results_38_TumortoT_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoT_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "T_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)


p1 <- plot_distance_boxplot(
  distance_results_38_TumortoT_Pre,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("A) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to T cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 150))  


p2 <- plot_distance_boxplot(
  distance_results_38_TumortoT_Post,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("B) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to T cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 150))  

# Radial layout plots

plot_radial_distance(
  distance_results_38_TumortoT_Pre,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to T cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


plot_radial_distance(
  distance_results_38_TumortoT_Post,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to T cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


# Analyzing the differences through t test
t.test(distance_results_38_TumortoT_Pre$T_cells, distance_results_38_TumortoT_Post$T_cells)


# Spatial network plots

visualize_spatial_network(
  posi_38,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "T_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "T_cells" = "#377EB8"),
  alpha = 0.7
) + ggtitle("Patient 38 - Spatial Network from Tumor to T cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),
      ncol = 2
    ))


visualize_spatial_network(
  posi_38_Pre,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "T_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "T_cells" = "#377EB8"),
  alpha = 0.7
)

visualize_spatial_network(
  posi_38_Post,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "T_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "T_cells" = "#377EB8"),
  alpha = 0.7
)


# Spatial Gradient Network IFN GAMMA

p3 <- visualize_spatial_gradient(
  spatial_data = posi_38_Pre,
  sample = "Patient_38",
  gradient_type = "T_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
) + scale_y_reverse() +
  ggtitle("C) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in T cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


p4 <- visualize_spatial_gradient(
  spatial_data = posi_38_Post,
  sample = "Patient_38",
  gradient_type = "T_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
)  + scale_y_reverse () + 
  ggtitle("D) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in T cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


# Correlation plots between distance Tumor to Myeloid and IFN GAMMA expression

merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoT_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)
merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoT_Post,
  by = "Newbarcode",
  all.x = FALSE 
)


result_correlation_38_Pre <- calculate_correlations(
  spatial_data = merged_38_IFN_Pre,
  distance_results = distance_results_38_TumortoT_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "T_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "E) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to T cells before Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Pre$estimate))
print(paste("P-value:", result_correlation_38_Pre$p_value))

p5 <- result_correlation_38_Pre$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 

result_correlation_38_Post <- calculate_correlations(
  spatial_data = merged_38_IFN_Post,
  distance_results = distance_results_38_TumortoT_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "T_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "F) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to T cells after Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Post$estimate))
print(paste("P-value:", result_correlation_38_Post$p_value))

p6 <- result_correlation_38_Post$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 



p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) 


pdf("Patient_38_STDistance_Tumor_T_cells.pdf", width = 28, height = 32)
p_all
dev.off()






########## B cells as targets ###############

# To calculate nearest distances to B cells

distance_results_38_TumortoB_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoB_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "B_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)


p1 <- plot_distance_boxplot(
  distance_results_38_TumortoB_Pre,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("A) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to B cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 350))  


p2 <- plot_distance_boxplot(
  distance_results_38_TumortoB_Post,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("B) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to B cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 350))  

# Radial layout plots

plot_radial_distance(
  distance_results_38_TumortoB_Pre,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to B cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


plot_radial_distance(
  distance_results_38_TumortoB_Post,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to B cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


# Analyzing the differences through t test
t.test(distance_results_38_TumortoB_Pre$B_cells, distance_results_38_TumortoB_Post$B_cells)



# Spatial network plots

visualize_spatial_network(
  posi_38,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "B_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "B_cells" = "#377EB8"),
  alpha = 0.7
) + ggtitle("Patient 38 - Spatial Network from Tumor to B cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),
      ncol = 2
    ))


visualize_spatial_network(
  posi_38_Pre,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "B_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "B_cells" = "#377EB8"),
  alpha = 0.7
)

visualize_spatial_network(
  posi_38_Post,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "B_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "B_cells" = "#377EB8"),
  alpha = 0.7
)


# Spatial Gradient Network IFN GAMMA

p3 <- visualize_spatial_gradient(
  spatial_data = posi_38_Pre,
  sample = "Patient_38",
  gradient_type = "B_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
) + scale_y_reverse() +
  ggtitle("C) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in B cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


p4 <- visualize_spatial_gradient(
  spatial_data = posi_38_Post,
  sample = "Patient_38",
  gradient_type = "B_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
)  +   scale_y_reverse() + 
  ggtitle("D) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in B cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


# Correlation plots between distance Tumor to B cells and IFN GAMMA expression

merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoB_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)
merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoB_Post,
  by = "Newbarcode",
  all.x = FALSE 
)


result_correlation_38_Pre <- calculate_correlations(
  spatial_data = posi_38_Pre,
  distance_results = distance_results_38_TumortoB_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "B_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "E) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to B cells before Denosumab"
)

print(paste("Correlation coefficient:", result_correlation_38_Pre$estimate))
print(paste("P-value:", result_correlation_38_Pre$p_value))

p5 <- result_correlation_38_Pre$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 

result_correlation_38_Post <- calculate_correlations(
  spatial_data = posi_38_Post,
  distance_results = distance_results_38_TumortoB_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "B_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "F) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to B cells after Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Post$estimate))
print(paste("P-value:", result_correlation_38_Post$p_value))

p6 <- result_correlation_38_Post$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 



p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) 


pdf("Patient_38_STDistance_Tumor_B_cells.pdf", width = 28, height = 32)
p_all
dev.off()





########## Endothelial cells as targets ###############

# To calculate nearest distances to Endothelial cells

distance_results_38_TumortoEndothelial_Pre <- calculate_nearest_distances(
  posi_38_Pre,
  reference_type = "tumor",   
  target_types = "Endothelial_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)

distance_results_38_TumortoEndothelial_Post <- calculate_nearest_distances(
  posi_38_Post,
  reference_type = "tumor",   
  target_types = "Endothelial_cells",
  x_col = "pxl_col_in_fullres",    
  y_col = "pxl_row_in_fullres",    
  id_col = "Newbarcode",
  type_col = "final_annotation"    
)


p1 <- plot_distance_boxplot(
  distance_results_38_TumortoEndothelial_Pre,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("A) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to Endothelial cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 125))  

p2 <- plot_distance_boxplot(
  distance_results_38_TumortoEndothelial_Post,
  id_col = "Newbarcode",
  show_points = TRUE,
  y_scale = "original",
  palette = "Dark2"
) + ggtitle("B) Patient 38 - Nearest Neighbor Distance Comparison \n from Tumor to Endothelial cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 125))  

# Radial layout plots

plot_radial_distance(
  distance_results_38_TumortoEndothelial_Pre,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to Endothelial cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))


plot_radial_distance(
  distance_results_38_TumortoEndothelial_Post,
  id_col = "Newbarcode",
  reference_type = "tumor",
  label_padding = 0.3,
  show_labels = TRUE,
  palette = "Dark2"
) + ggtitle("Patient 38 - Radial Layout from Tumor to Endothelial cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))



# Analyzing the differences through t test
t.test(distance_results_38_TumortoEndothelial_Pre$Endothelial_cells, distance_results_38_TumortoEndothelial_Post$Endothelial_cells)



# Spatial network plots

visualize_spatial_network(
  posi_38,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Endothelial_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Endothelial_cells" = "#377EB8"),
  alpha = 0.7
) + ggtitle("Patient 38 - Spatial Network from Tumor to Endothelial cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5),
      ncol = 2
    ))


visualize_spatial_network(
  posi_38_Pre,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Endothelial_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Endothelial_cells" = "#377EB8"),
  alpha = 0.7
)

visualize_spatial_network(
  posi_38_Post,
  sample = "Patient_38",
  reference_type = "tumor",
  target_type = "Endothelial_cells",
  x_col = "pxl_col_in_fullres",
  y_col = "pxl_row_in_fullres",
  type_col = "final_annotation",
  color_palette = c("tumor" = "#90ee90", "Endothelial_cells" = "#377EB8"),
  alpha = 0.7
)


# Spatial Gradient Network IFN GAMMA

p3 <- visualize_spatial_gradient(
  spatial_data = posi_38_Pre,
  sample = "Patient_38",
  gradient_type = "Endothelial_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
) +   scale_y_reverse() + 
  ggtitle("C) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in Endothelial cells before Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


p4 <- visualize_spatial_gradient(
  spatial_data = posi_38_Post,
  sample = "Patient_38",
  gradient_type = "Endothelial_cells",
  fixed_type = "tumor",
  expression_col = "ONLY_EXP_UP_IFN_GAMMA",
  type_col = "final_annotation",
  fixed_color = "#CCCCCC",
  line_color = "#444444",
  gradient_palette = "C", 
  point_size = 1.5,
  point_alpha = 0.9
)  +   scale_y_reverse() + 
  ggtitle("D) Patient 38 - Spatial Gradient Network - Tumor as reference and INTERFERON GAMMA (Bulk-RNAseq) \n expression in Endothelial cells after Denosumab") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  ) 


# Correlation plots between distance Tumor to Endothelial cells and IFN GAMMA expression

merged_38_IFN_Pre <- merge(
  posi_38_Pre,
  distance_results_38_TumortoEndothelial_Pre,
  by = "Newbarcode",
  all.x = FALSE 
)
merged_38_IFN_Post <- merge(
  posi_38_Post,
  distance_results_38_TumortoEndothelial_Post,
  by = "Newbarcode",
  all.x = FALSE 
)


result_correlation_38_Pre <- calculate_correlations(
  spatial_data = posi_38_Pre,
  distance_results = distance_results_38_TumortoEndothelial_Pre,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Endothelial_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "E) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to Endothelial cells before Denosumab"
)

print(paste("Correlation coefficient:", result_correlation_38_Pre$estimate))
print(paste("P-value:", result_correlation_38_Pre$p_value))

p5 <- result_correlation_38_Pre$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 

result_correlation_38_Post <- calculate_correlations(
  spatial_data = posi_38_Post,
  distance_results = distance_results_38_TumortoEndothelial_Post,
  spatial_feature = "ONLY_EXP_UP_IFN_GAMMA",
  distance_metric = "Endothelial_cells",
  method = "pearson",
  plot = TRUE,
  plot_title = "F) Patient 38 - Correlation between INTERFERON GAMMA (Bulk-RNAseq) Expression \n and distance from Tumor to Endothelial cells after Denosumab"
) 

print(paste("Correlation coefficient:", result_correlation_38_Post$estimate))
print(paste("P-value:", result_correlation_38_Post$p_value))

p6 <- result_correlation_38_Post$plot + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  legend.title = element_text(size = 16, face = "bold"),
  axis.text.x  = element_text(size = 14),  
  axis.text.y  = element_text(size = 14),  
  axis.title.x = element_text(size = 16),  
  axis.title.y = element_text(size = 16)
) 



p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) 


pdf("Patient_38_STDistance_Tumor_Endothelial_cells.pdf", width = 28, height = 32)
p_all
dev.off()

