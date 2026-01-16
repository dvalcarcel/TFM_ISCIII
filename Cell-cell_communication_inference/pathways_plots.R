library(ggplot2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(stringr)
library(readxl)
library(ggsci)
library(forcats)


setwd("PATH")

df <- read_excel("All_Patients_CellChat_Interactions_50micras.xls")
df <- as.data.frame(df)

#Adding pre and post conditions
df_mod <- df |> 
  mutate(Condition = ifelse(grepl("Pre", source), "Pre", "Post"))

#Remove endings from cell identities since there is one column for patient and another for condition
df_mod <- df_mod %>%
  mutate(source = str_remove(source, "(_Pre|_Post|_Pre_ptx4|_Post_ptx4|_Pre_ptx26|_Post_ptx26|_Pre_ptx41|_Post_ptx41|_Pre_ptx67|_Post_ptx67)$"))
df_mod <- df_mod %>%
  mutate(target = str_remove(target, "(_Pre|_Post|_Pre_ptx4|_Post_ptx4|_Pre_ptx26|_Post_ptx26|_Pre_ptx41|_Post_ptx41|_Pre_ptx67|_Post_ptx67)$"))


#Filtering specific cell idents
cell_idents <- c(
  "Stromal",
  "Endothelial_cells",
  "Myeloid",
  "tumor",
  "T_cells",
  "B_cells"
)

df_mod <- df_mod[df_mod$source %in% cell_idents &df_mod$target %in% cell_idents,]

# Calculate TotalProb column to sum probabilities and transform to plot
df_long_adj <- df_mod %>%
  mutate(prob = as.numeric(prob),
         CellType = target) %>%
  group_by(Patient, source, CellType, Condition, interaction_name, pathway_name) %>%
  summarise(TotalProb = sum(prob, na.rm = TRUE), .groups = "drop")

write.table(df_long_adj, file="All_Patients_CellChat_TotalProb_by_Pathway_50um_filtered.txt", sep = '\t')


#### Defining the color for each pathway #####
pathway_colors <- c(
  APP="#e41a1c", COLLAGEN="#ff7f0e", LAMININ="#4daf4a", PECAM1="#17becf",
  SPP1="#9467bd", FN1="#8c564b", GAP="#f781bf", THBS="#999999",
  Netrin="#66c2a5", CDH="#fc8d62", MIF="#8da0cb", SEMA3="#e78ac3",
  COMPLEMENT="#a6d854", CXCL="#ffd92f", ADIPONECTIN="#e5c494",
  CCL="#b3b3b3", IL16="#1b9e77", SELL="#d95f02", GRN="#a9961d",
  PTPR="#e7298a", PECAM2="#66a61e", TENASCIN="#e6ab02", ADGRG="#a6761d",
  CD96="#666666", CD99="#8dd3c7", IGFBP="#ffffb3", CDH5="#bebada",
  ESAM="#fb8072", NOTCH="#80b1d3", ANGPTL="#fdb462", CD46="#b3de69",
  EPHB="#fccde5", GAS="#d9d9d9", JAM="#bc80bd", MK="#ccebc5",
  PDGF="#ffed6f", VEGF="#1f78b4", PERIOSTIN="#33a02c", PTPRM="#fb9a99",
  ADGRA="#a6cee3", SEMA4="#b2df8a", SEMA6="#fdbf6f", TGFb="#cab2d6",
  ANNEXIN="#ffff99", PLAU="#6a3d9a", THY1="#ff7f7f", VCAM="#b15928",
  MMP="#e31a1c", VWF="#6a51a3", HSPG="#2ca25f", AGRN="#e7298a",
  EPHA="#3182bd", ACTIVIN="#31a354", ADGRB="#756bb1", RELN="#636363",
  CSF="#43a2ca", PARs="#fb6a4a", RA="#238b45", TWEAK="#a1d99b",
  Prostaglandin="#bcbd22", ADGRE="#9467bd", PTN="#17becf",
  DESMOSOME="#9edae5", CysLTs="#8c564b", MPZ="#c49c94", CypA="#e377c2",
  SEMA5="#7f7f7f", CDH1="#aec7e8", Desmosterol="#ffbb78", KIT="#98df8a",
  ICAM="#ff9896", NECTIN="#c5b0d5", CEACAM="#c7c7c7", OCLN="#dbdb8d",
  WNT="#9edae5", PCDH="#17becf", GALECTIN="#1f77b4", BAFF="#ff7f0e",
  CD6="#2ca02c", CD39="#fb8072", CALCR="#9467bd", APRIL="#8c564b",
  CADM="#e377c2", BMP="#7f7f7f", CLDN="#bcbd22", Cholesterol="#79becf",
  DHT="#aec7e8", EGF="#ffbb78", VISFATIN="#98df8a", ncWNT="#ff9896",
  CSPG4="#c5b0d5"
)




################# Patient 38 ###############################

df_38 <- df_long_adj[df_long_adj$Patient==38, ]

###### Myeloid ######

df_38_myel <- df_38[df_38$source %in% "Myeloid", ]

y_max_myel <- 0.0015

png("Patient_38_Myeloid_cells_50µm_Pre_vs_Post.png", width = 9000, height = 4500, res = 800)
df_38_myel %>%
  mutate(Condition = factor(Condition, levels = c("Pre", "Post"))) %>%
  ggplot(aes(x = reorder(CellType, TotalProb, FUN = sum), 
             y = TotalProb, fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = pathway_colors) +  
  labs(
    x = "Target cell type",
    y = "Sum of interaction probabilities",
    title = "D) Patient 38 - Myeloid cell-cell interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  coord_cartesian(ylim = c(0, y_max_myel)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,  face = "bold")
                              
  )
dev.off()



###### T cells ######

df_38_t <- df_38[df_38$source %in% "T_cells", ]

y_max_t <- 0.01

p1 <- df_38_t %>%
  mutate(Condition = factor(Condition, levels = c("Pre", "Post"))) %>%
  ggplot(aes(x = reorder(CellType, TotalProb, FUN = sum), 
             y = TotalProb, fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = pathway_colors) +  
  labs(
    x = "Target cell type",
    y = "Sum of interaction probabilities",
    title = "Patient 38 - T cell-cell interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  coord_cartesian(ylim = c(0, y_max_t)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,  face = "bold")
    
  )


###### tumor ######

df_38_tumor <- df_38[df_38$source %in% "tumor", ]

y_max_tumor <- 0.01


p2 <- df_38_tumor %>%
  mutate(Condition = factor(Condition, levels = c("Pre", "Post"))) %>%
  ggplot(aes(x = reorder(CellType, TotalProb, FUN = sum), 
             y = TotalProb, fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = pathway_colors) +  
  labs(
    x = "Target cell type",
    y = "Sum of interaction probabilities",
    title = "Patient 38 - Tumor cell-cell interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  coord_cartesian(ylim = c(0, y_max_t)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,  face = "bold")
    
  )


###### Stromal ######

df_38_stromal <- df_38[df_38$source %in% "Stromal", ]

y_max_stromal <- 0.01

p3 <- df_38_stromal %>%
  mutate(Condition = factor(Condition, levels = c("Pre", "Post"))) %>%
  ggplot(aes(x = reorder(CellType, TotalProb, FUN = sum), 
             y = TotalProb, fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = pathway_colors) +  
  labs(
    x = "Target cell type",
    y = "Sum of interaction probabilities",
    title = "Patient 38 - Stromal cell-cell interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  coord_cartesian(ylim = c(0, y_max_stromal)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,  face = "bold")
    
  )


###### Endothelial ######

df_38_endothelial <- df_38[df_38$source %in% "Endothelial_cells", ]

y_max_endothelial <- 0.01

p4 <- df_38_endothelial %>%
  mutate(Condition = factor(Condition, levels = c("Pre", "Post"))) %>%
  ggplot(aes(x = reorder(CellType, TotalProb, FUN = sum), 
             y = TotalProb, fill = pathway_name)) +
  geom_col(position = "stack") +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = pathway_colors) +  
  labs(
    x = "Target cell type",
    y = "Sum of interaction probabilities",
    title = "Patient 38 - Endothelial cell-cell interactions by pathway (50µm)",
    fill = "Pathway"
  ) +
  coord_cartesian(ylim = c(0, y_max_endothelial)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16,  face = "bold")
    
  )


pdf("Patient_38_T_tumor_stromal_endothelial__50µm_Pre_vs_Post.pdf", height = 15, width = 25)
p1+p2+p3+p4
dev.off()
