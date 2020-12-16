# 152 libraries across 102 experiments
# ENCODE 3 (34) and ENCODE 4 (68) experiments including tissue, differentiated cells, and cell lines

library(ggplot2)
library(tidyr)
library(plyr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(Seurat)

setwd("/Users/liz/Documents/miRNA/ee282")

#################### PCA: log2((CPM > 2)+1) #################### 
load("rdata/cpm_2_geneid.rda")
load("rdata/cpm_2_geneid_ti.rda")
load("rdata/cpm_2_geneid_cl_dcl.rda")

# Samples are in rows for PCA, transpose data frame
pca_func <- function(cpm_2) {
  cpm_log2 <- t(log2(cpm_2+1))
  cpm_log2_pca <- prcomp(cpm_log2)
  cpm_log2_pca_out <- as.data.frame(cpm_log2_pca$x)
  cpm_log2_pca_out$var_explained = (cpm_log2_pca$sdev/sum(cpm_log2_pca$sdev))*100
  cpm_log2_pca_out$miRNA_ENCODE_Data_ID = rownames(cpm_log2_pca_out)
  return(cpm_log2_pca_out)
}

cpm_2_pca = pca_func(cpm_2_geneid)
cpm_2_pca_ti = pca_func(cpm_2_geneid_ti)
cpm_2_pca_cl_dcl = pca_func(cpm_2_geneid_cl_dcl)

# Merge metadata with matrix by miRNA_ENCODE_Data_ID
cpm_2_pca_plus_metadata = join(cpm_2_pca, metadata)
cpm_2_pca_ti_plus_metadata = join(cpm_2_pca_ti, metadata_ti)
cpm_2_pca_cl_dcl_plus_metadata = join(cpm_2_pca_cl_dcl, metadata_cl_dcl)

#################### UMAP on all samples  #################### 
load("rdata/counts_0_geneIDs.rda")
table(colnames(counts_0_geneid) == metadata$miRNA_ENCODE_Data_ID) 

mirna <- CreateSeuratObject(counts = counts_0_geneid, project = "miRNA", min.cells = 0, min.features = 0)
mirna
mirna@meta.data$BroaderTissueType <- metadata$BroaderTissueType
mirna@meta.data$TissueType <- metadata$TissueType
mirna@meta.data$SampleType <- metadata$SampleType
mirna <- NormalizeData(mirna)
mirna <- FindVariableFeatures(mirna, selection.method = "vst", nfeatures = 100)
all.genes <- rownames(mirna)
mirna <- ScaleData(mirna, features = all.genes)
mirna <- RunPCA(mirna, features = VariableFeatures(object = mirna), npcs = 99)
ElbowPlot(mirna, ndims  = 99)
mirna <- RunUMAP(mirna, dims = 1:10)
nb.cols <- length(unique(mirna$BroaderTissueType))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)


#################### Plot PCA: color by sample #################### 
pca_plot_sample <- function(cpm_2_pca_and_metadata) {
  percentage <- round(cpm_2_pca_and_metadata$var_explained,2)
  percentage_formatted <- paste0("PC", paste(as.character(1:length(cpm_2_pca_and_metadata$var_explained))), " (", paste(as.character(percentage)), "%)")
  nb.cols2 <- length(unique(cpm_2_pca_and_metadata$SampleType))
  mycolors2 <- colorRampPalette(brewer.pal(3, "Set1"))(nb.cols2)
  p_sample <- ggplot(cpm_2_pca_and_metadata,aes(x=PC1,y=PC2,color=SampleType)) + geom_point(size = 1) + theme_bw() + 
    xlab(percentage_formatted[1]) + ylab(percentage_formatted[2]) +
    ggtitle("PCA of 152 microRNA-seq libraries\nlabeled by sample type")+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8, face="bold"),
          legend.title=element_blank(), 
          legend.position = "none",
         # legend.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  + scale_color_manual(values=mycolors2)
  return(p_sample)
  
}

#################### Compare PCA and UMAP #################### 
p1 = pca_plot_sample(cpm_2_pca_plus_metadata)

p2 = DimPlot(mirna, reduction = "umap", pt.size = 1, 
        group.by = "SampleType",
        label=F, cols = c("#e41a19","#357eb7","#4cae48")) + 
  ggtitle("UMAP of 152 microRNA-seq libraries\nlabeled by sample type")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  

png(file = "figures/pca_umap_allSamples_SampleType.png",
    width = 7, 
    height = 4, units = "in", res = 1000)
ggarrange(p1,p2,common.legend = TRUE,legend="bottom")
dev.off()

#################### PCA: color by tissue type #################### 
pca_plot_broaderTissue <- function(cpm_2_pca_and_metadata) {
  # Get percent variance explained for plotting
  percentage <- round(cpm_2_pca_and_metadata$var_explained,2)
  percentage_formatted <- paste0("PC", paste(as.character(1:length(cpm_2_pca_and_metadata$var_explained))), " (", paste(as.character(percentage)), "%)")
  # Define the number of colors you want
  nb.cols3 <- length(unique(cpm_2_pca_and_metadata$BroaderTissueType))
  mycolors3 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols3)
  
  p_broadTissuetype <- ggplot(cpm_2_pca_and_metadata,aes(x=PC1,y=PC2,color=BroaderTissueType, label=BroaderTissueType)) +
    geom_point(size = 1) +
    geom_text_repel(data=subset(cpm_2_pca_and_metadata, 
                                grepl("Heart_LeftAtrium|Myocyte", cpm_2_pca_and_metadata$TissueType)), 
                                 size = 3,box.padding = 0.35, 
                                 point.padding = 0.5,
                                 segment.color = 'grey50', show.legend=FALSE) + theme_bw() + 
    xlab(percentage_formatted[1]) + ylab(percentage_formatted[2]) +
    ggtitle("PCA of 152 microRNA-seq libraries\nlabeled by derived/simulated tissue type")+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8, face="bold"),
          legend.title=element_text(size=8), 
          legend.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  + scale_color_manual(values=mycolors3)
  return(p_broadTissuetype)
}

#################### Compare PCA and UMAP #################### 
p3 = pca_plot_broaderTissue(cpm_2_pca_plus_metadata)

p4 = DimPlot(mirna, reduction = "umap", pt.size = 1, group.by = "BroaderTissueType", 
  label=T, label.size = 2, repel = T, cols = mycolors)+
  ggtitle("UMAP of 152 microRNA-seq libraries\nlabeled by derived/simulated tissue type")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

png(file = "figures/pca_umap_allSamples_TissueType.png",
    width = 7, 
    height = 5.5, units = "in", res = 1000)
ggarrange(p3,p4,common.legend = TRUE,legend="bottom")
dev.off()


################ UMAP and PCA - tissues only ################ 
pca_plot_broaderTissue2 <- function(cpm_2_pca_and_metadata) {
  # Get percent variance explained for plotting
  percentage <- round(cpm_2_pca_and_metadata$var_explained,2)
  percentage_formatted <- paste0("PC", paste(as.character(1:length(cpm_2_pca_and_metadata$var_explained))), " (", paste(as.character(percentage)), "%)")
  # Define the number of colors you want
  nb.cols3 <- length(unique(cpm_2_pca_and_metadata$BroaderTissueType))
  mycolors3 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols3)
  
  p_broadTissuetype <- ggplot(cpm_2_pca_and_metadata,aes(x=PC1,y=PC2,color=BroaderTissueType, label=BroaderTissueType)) +
    geom_point(size = 1) + theme_bw() + 
    xlab(percentage_formatted[1]) + ylab(percentage_formatted[2]) +
    ggtitle("PCA of 86 tissue microRNA-seq libraries\nlabeled by derived/simulated tissue type")+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title = element_text(size = 8, face="bold"),
          legend.title=element_text(size=8), 
          legend.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  + scale_color_manual(values=mycolors3)
  return(p_broadTissuetype)
}

p5 = pca_plot_broaderTissue2(cpm_2_pca_ti_plus_metadata)

load("rdata/counts_0_geneids_tissues.rda")
table(colnames(counts_0_geneid_ti) == metadata_ti$miRNA_ENCODE_Data_ID) # sanity check

mirna <- CreateSeuratObject(counts = counts_0_geneid_ti, project = "miRNA", min.cells = 0, min.features = 0)
mirna
mirna@meta.data$BroaderTissueType <- metadata_ti$BroaderTissueType
mirna@meta.data$TissueType <- metadata_ti$TissueType
mirna@meta.data$SampleType <- metadata_ti$SampleType
mirna <- NormalizeData(mirna)
mirna <- FindVariableFeatures(mirna, selection.method = "vst", nfeatures = 100)
all.genes <- rownames(mirna)
mirna <- ScaleData(mirna, features = all.genes)
mirna <- RunPCA(mirna, features = VariableFeatures(object = mirna), npcs = 85)
ElbowPlot(mirna, ndims  = 85)
mirna <- RunUMAP(mirna, dims = 1:5)
nb.cols4 <- length(unique(mirna$BroaderTissueType))
mycolors4 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols4)


p6 = DimPlot(mirna, reduction = "umap", pt.size = 1, repel = T,
             group.by = "BroaderTissueType", label=T, label.size = 2, cols = mycolors4)+ 
  ggtitle("UMAP of 86 tissue microRNA-seq libraries\nlabeled by derived/simulated tissue type")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

png(file = "figures/pca_umap_TissueType_Tissues.png",
    width = 7, 
    height = 4.5, units = "in", res = 1000)
ggarrange(p5,p6,common.legend = TRUE,legend="bottom")
dev.off()

################ UMAP and PCA - cell lines/differentiated cell lines only ################ 
p7 = pca_plot_broaderTissue2(cpm_2_pca_cl_dcl_plus_metadata) 
p7 = p7 + ggtitle("PCA of 66 cell line/differentiated cell line\nmicroRNA-seq libraries labeled by derived/simulated tissue type")

load("rdata/counts_0_geneids_cls_dcls.rda")
table(colnames(counts_0_geneid_cl_dcl) == metadata_cl_dcl$miRNA_ENCODE_Data_ID) # sanity check

mirna <- CreateSeuratObject(counts = counts_0_geneid_cl_dcl, project = "miRNA", min.cells = 0, min.features = 0)
mirna
mirna@meta.data$BroaderTissueType <- metadata_cl_dcl$BroaderTissueType
mirna@meta.data$TissueType <- metadata_cl_dcl$TissueType
mirna@meta.data$SampleType <- metadata_cl_dcl$SampleType
mirna <- NormalizeData(mirna)
mirna <- FindVariableFeatures(mirna, selection.method = "vst", nfeatures = 100)
all.genes <- rownames(mirna)
mirna <- ScaleData(mirna, features = all.genes)
mirna <- RunPCA(mirna, features = VariableFeatures(object = mirna), npcs = 65)
ElbowPlot(mirna, ndims  = 65)
mirna <- RunUMAP(mirna, dims = 1:12)
nb.cols4 <- length(unique(mirna$BroaderTissueType))
mycolors4 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols4)


p8 = DimPlot(mirna, reduction = "umap", pt.size = 1, repel = T,
             group.by = "BroaderTissueType", label=T, label.size = 2, cols = mycolors4)+ 
  ggtitle("UMAP of 66 cell line/differentiated cell line\nmicroRNA-seq libraries labeled by derived/simulated tissue type")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 8),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

png(file = "figures/pca_umap_TissueType_CLsDCLs.png",
    width = 7, 
    height = 4.5, units = "in", res = 1000)
ggarrange(p7,p8,common.legend = TRUE,legend="bottom")
dev.off()

########
########

