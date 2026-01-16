library(Seurat)
library(scProportionTest)
library(ggplot2)

setwd("PATH")

load("Patient_4_26_RCTD_snRNAseq.rda")

Patient_4_26_RCTD$Patient_Denosumab <- paste(Patient_4_26_RCTD$Patient,Patient_4_26_RCTD$Denosumab, sep = "_")

prop_test <- sc_utils(Patient_4_26_RCTD)

prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "ptx4_Pre", sample_2 = "ptx4_Post",
  sample_identity = "Patient_Denosumab",
  n_permutations = 1000
)

p1 <- permutation_plot(prop_test) +
  ggtitle("A) Patient 4 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )



load("Patient_06_RCTD_snRNAseq.rda")

prop_test <- sc_utils(Patient_06_RCTD)

prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "Pre", sample_2 = "Post",
  sample_identity = "Denosumab",
  n_permutations = 1000
)

p2 <- permutation_plot(prop_test) +
  ggtitle("B) Patient 6 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )

png("Patient_4_6_scProportionTest.png", width = 12000, height = 4000, res = 650)
p1+p2
dev.off()



load("Patient_38_RCTD_snRNAseq.rda")

prop_test <- sc_utils(Patient_38_RCTD)

prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "Pre", sample_2 = "Post",
  sample_identity = "Denosumab",
  n_permutations = 1000
)

png("Patient_38_scProportionTest.png", width = 6000, height = 4000, res = 650)
permutation_plot(prop_test) +
  ggtitle("Patient 38 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )
dev.off()



load("Patient_41_67_RCTD_snRNAseq.rda")

Patient_41_67_RCTD$Patient_Denosumab <- paste(Patient_41_67_RCTD$Patient,Patient_41_67_RCTD$Denosumab, sep = "_")

prop_test <- sc_utils(Patient_41_67_RCTD)


prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "ptx41_Pre", sample_2 = "ptx41_Post",
  sample_identity = "Patient_Denosumab",
  n_permutations = 1000
)

p3 <- permutation_plot(prop_test) +
  ggtitle("C) Patient 41 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )

prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "ptx67_Pre", sample_2 = "ptx67_Post",
  sample_identity = "Patient_Denosumab",
  n_permutations = 1000
)

p4 <- permutation_plot(prop_test) +
  ggtitle("D) Patient 67 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )

png("Patient_41_67_scProportionTest.png", width = 12000, height = 4000, res = 650)
p3+p4
dev.off()


load("Patient_61_RCTD_snRNAseq.rda")

prop_test <- sc_utils(Patient_61_RCTD)

prop_test <- permutation_test(
  prop_test, cluster_identity = "final_annotation",
  sample_1 = "Pre", sample_2 = "Post",
  sample_identity = "Denosumab",
  n_permutations = 1000
)

p5 <- permutation_plot(prop_test) +
  ggtitle("E) Patient 61 – scProportionTest") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14),
    axis.text.x  = element_text(size = 14),  
    axis.text.y  = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16)
  )


png("Patient_61_scProportionTest.png", width = 6000, height = 4000, res = 650)
p5
dev.off()
