setwd("PATH")
library("DESeq2")
library("stats")
library("ggplot2")
library("ggrepel")
library("dplyr")

#CountTable derived from FindAllMarkers pb matrix
countTable = read.table("B_cell_pb_matrix.csv",header=TRUE, row.names=1, check.names = FALSE,sep = ',')
designMatrix = read.table('B_cell_metadata.csv', header=TRUE, row.names=1,sep = ',')
libType = c("oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor","oneFactor")
condition = factor(designMatrix$Condition)
patient = factor(designMatrix$Patient)
condition <- relevel(condition, "Pre")
experiment_design=data.frame(row.names = colnames(countTable),condition,patient,libType)

colnames(countTable) == rownames(experiment_design)

cds <- DESeqDataSetFromMatrix(countData =countTable, colData=experiment_design, design=~patient + condition)
cds_DESeqED <- DESeq(cds,parallel = TRUE)
res <- results(cds_DESeqED,parallel = TRUE,alpha = 0.05, pAdjustMethod = "BH")
write.table(res,file = "Bcell_Post_vs_Pre.differentialExpression.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
resSig <- subset(res, padj < 0.05)
write.table(resSig,file = "Bcell_Post_vs_Pre.sigDEGs.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")
normalizedReadCounts = counts(cds_DESeqED,normalized=TRUE)
write.table(normalizedReadCounts,file = "Bcell_Post_vs_Pre.normalizedCounts.txt",row.names = TRUE,col.names = NA,append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".")


## Volcano plot with labels

Bcell_Post_vs_Pre <- read.table("/Users/dvalcarcel/Documents/TFM/Analysis_16um/RCTD_deconvolution/RCTC_snRNAseq/v2/PB_final_annotation/Bcell_Post_vs_Pre.differentialExpression.txt", header=TRUE, sep = '\t')
Bcell_Post_vs_Pre <- Bcell_Post_vs_Pre[complete.cases(Bcell_Post_vs_Pre), ]

Bcell_Post_vs_Pre$Type <- "NONE"
Bcell_Post_vs_Pre$Type[Bcell_Post_vs_Pre$padj < 0.05 & Bcell_Post_vs_Pre$log2FoldChange < 0] <- "DOWN"
Bcell_Post_vs_Pre$Type[Bcell_Post_vs_Pre$padj < 0.05 & Bcell_Post_vs_Pre$log2FoldChange > 0] <- "UP"
cols <- c("NONE" = "#474657", "DOWN" = "#157ded", "UP" = "#f2800f")

Bcell_Post_vs_Pre$genelabels <- ""

# To make a column with labels according to the threshold
Bcell_Post_vs_Pre$label <- ifelse(
  Bcell_Post_vs_Pre$padj < 0.05,
  Bcell_Post_vs_Pre$X,
  NA
)

png("DESeq2_B_Post_vs_Pre.png", width = 10000, height = 8000, res = 650)

p2 <- ggplot(Bcell_Post_vs_Pre, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 0.75, na.rm = TRUE) +
  geom_text_repel(aes(label = label), size = 8, max.overlaps = 30, na.rm = TRUE) +
  theme_bw(base_size = 18) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20)) +
  xlab(expression(log[2]("Fold Change"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = 1.3, colour = "red4", linetype = "dashed") + 
  geom_vline(xintercept = 0.0, colour = "#474657", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  ggtitle("C) DESeq2 Volcano Plot - B cells Pre vs Post") 

dev.off()
