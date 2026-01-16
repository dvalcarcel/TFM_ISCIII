library(loupeR)
library(stringr)
library(Seurat)

setwd("PATH")

# We need to import the annotation of Denosumab clustes (i.e. pre1, pre2 and post)
# from LoupeBrowser to the Seurat object in order to perform semla analysis

#Function to transform 8um spots into 16um.

TransformBarcodes<-function(Barcodes,SizeOriginal,SizeNew)
{
  SizeO<-str_pad(paste0(SizeOriginal,"um"),5,pad="0")
  SizeN<-str_pad(paste0(SizeNew,"um"),5,pad="0")
  
  Barcodes<-strsplit(Barcodes,"[_-]")
  Barcodes<-data.frame(s="s",size=SizeO,
                       X=as.numeric(sapply(Barcodes,function(X){return(X[3])})),
                       Y=as.numeric(sapply(Barcodes,function(X){return(X[4])})),
                       end="-1")
  
  
  nBinsSide <- SizeOriginal / SizeNew
  
  Xmins <- Barcodes$X * nBinsSide
  Xmaxs <- (Barcodes$X * nBinsSide) + (nBinsSide - 1)
  Ymins <- Barcodes$Y * nBinsSide
  Ymaxs <- (Barcodes$Y * nBinsSide) + (nBinsSide - 1)
  
  Result<-lapply(1:nrow(Barcodes),function(jj)
  {
    Rx <- expand.grid(x = Xmins[jj]:Xmaxs[jj], y = Ymins[jj]:Ymaxs[jj])
    Rx$X <- Barcodes$X[jj]
    Rx$Y <- Barcodes$Y[jj]
    return(Rx)
  })
  
  Result <- do.call(rbind, Result)
  
  OldBC<-paste0("s_",SizeO,"_",str_pad(Result$X,5,pad="0"),"_",str_pad(Result$Y,5,pad="0"),"-1")
  NewBC<-paste0("s_",SizeN,"_",str_pad(Result$x,5,pad="0"),"_",str_pad(Result$y,5,pad="0"),"-1")
  
  Result<-data.frame(Original=OldBC,Transformed=NewBC)
  #Result<-split(Result$Transformed,Result$Original)
  
  return(Result)
  
}

### Patient_38

#Importing the csv annotation from LoupeBrowser
Patient_38_metadata <- read.csv("Denosumab_cluster_38.csv")
Patient_38_8um_bc <- Patient_38_metadata$Barcode
Patient_38_16um_bc <- TransformBarcodes(Patient_38_8um_bc,SizeOriginal = 8, SizeNew = 16)
rownames(Patient_38_metadata) <- Patient_38_16um_bc$Transformed


load("Patient_38_RCTD_snRNAseq.rda")


Patient_38_metadata_16um <- Patient_38_metadata[colnames(Patient_38_RCTD), ]
Patient_38_RCTD <- AddMetaData(object = Patient_38_RCTD, metadata = Patient_38_metadata_16um$Denosumab_cluster, col.name = "Denosumab_cluster")

write.csv(Patient_38_metadata_16um, file = "Patient_38_16um_metadata.csv")
save(Patient_38_RCTD, file="Patient_38_RCTD_snRNAseq.rda")

pdf("Patient_38_Denosumab_cluster.pdf")
SpatialDimPlot(Patient_38_RCTD, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
dev.off()


### Patient_61


Patient_61_metadata <- read.csv("Denosumab_cluster_61.csv")
Patient_61_8um_bc <- Patient_61_metadata$Barcode
Patient_61_16um_bc <- TransformBarcodes(Patient_61_8um_bc,SizeOriginal = 8, SizeNew = 16)
rownames(Patient_61_metadata) <- Patient_61_16um_bc$Transformed


load("Patient_61_RCTD_snRNAseq.rda")


Patient_61_metadata_16um <- Patient_61_metadata[colnames(Patient_61_RCTD), ]
Patient_61_RCTD <- AddMetaData(object = Patient_61_RCTD, metadata = Patient_61_metadata_16um$Denosumab_cluster, col.name = "Denosumab_cluster")

write.csv(Patient_61_metadata_16um, file = "Patient_61_16um_metadata.csv")
save(Patient_61_RCTD, file="Patient_61_RCTD_snRNAseq.rda")

pdf("Patient_61_Denosumab_cluster.pdf")
SpatialDimPlot(Patient_61_RCTD, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
dev.off()




### Patient_06

Patient_06_metadata <- read.csv("Denosumab_cluster_06.csv")
Patient_06_8um_bc <- Patient_06_metadata$Barcode
Patient_06_16um_bc <- TransformBarcodes(Patient_06_8um_bc,SizeOriginal = 8, SizeNew = 16)
rownames(Patient_06_metadata) <- Patient_06_16um_bc$Transformed


load("Patient_06_RCTD_snRNAseq.rda")


Patient_06_metadata_16um <- Patient_06_metadata[colnames(Patient_06_RCTD), ]
Patient_06_RCTD <- AddMetaData(object = Patient_06_RCTD, metadata = Patient_06_metadata_16um$Denosumab_cluster, col.name = "Denosumab_cluster")

write.csv(Patient_06_metadata_16um, file = "Patient_06_16um_metadata.csv")
save(Patient_06_RCTD, file="Patient_06_RCTD_snRNAseq.rda")

pdf("Patient_06_Denosumab_cluster.pdf")
SpatialDimPlot(Patient_06_RCTD, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
dev.off()

### Patient_41_67

Patient_41_67_metadata <- read.csv("Denosumab_cluster_41_67.csv")
Patient_41_67_8um_bc <- Patient_41_67_metadata$Barcode
Patient_41_67_16um_bc <- TransformBarcodes(Patient_41_67_8um_bc,SizeOriginal = 8, SizeNew = 16)
rownames(Patient_41_67_metadata) <- Patient_41_67_16um_bc$Transformed


load("Patient_41_67_RCTD_snRNAseq.rda")


Patient_41_67_metadata_16um <- Patient_41_67_metadata[colnames(Patient_41_67_RCTD), ]
Patient_41_67_RCTD <- AddMetaData(object = Patient_41_67_RCTD, metadata = Patient_41_67_metadata_16um$Denosumab_cluster, col.name = "Denosumab_cluster")

write.csv(Patient_41_67_metadata_16um, file = "Patient_41_67_16um_metadata.csv")
save(Patient_41_67_RCTD, file="Patient_41_67_RCTD_snRNAseq.rda")


pdf("Patient_41_67_Denosumab_cluster.pdf")
SpatialDimPlot(Patient_41_67_RCTD, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
dev.off()




### Patient_4_26

Patient_4_26_metadata <- read.csv("Denosumab_cluster_4_26.csv")
Patient_4_26_8um_bc <- Patient_4_26_metadata$Barcode
Patient_4_26_16um_bc <- TransformBarcodes(Patient_4_26_8um_bc,SizeOriginal = 8, SizeNew = 16)
rownames(Patient_4_26_metadata) <- Patient_4_26_16um_bc$Transformed


load("Patient_4_26_RCTD_snRNAseq.rda")


Patient_4_26_metadata_16um <- Patient_4_26_metadata[colnames(Patient_4_26_RCTD), ]
Patient_4_26_RCTD <- AddMetaData(object = Patient_4_26_RCTD, metadata = Patient_4_26_metadata_16um$Denosumab_cluster, col.name = "Denosumab_cluster")

write.csv(Patient_4_26_metadata_16um, file = "Patient_4_26_16um_metadata.csv")
save(Patient_4_26_RCTD, file="Patient_4_26_RCTD_snRNAseq.rda")

pdf("Patient_4_26_Denosumab_cluster.pdf")
SpatialDimPlot(Patient_4_26_RCTD, group.by ="Denosumab_cluster", image.alpha=0, pt.size.factor=3.2)
dev.off()

