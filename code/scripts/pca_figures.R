# 152 libraries across 102 experiments
# ENCODE 3 (34) and ENCODE 4 (68) experiments including tissue, differentiated cells, and cell lines

library(ggplot2)
library(tidyr)
library(plyr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("/Users/liz/Documents/miRNA/ee282")

#################### Load metadata and data #################### s
load("rdata/cpm_2_geneid.rda")
load("rdata/cpm_2_geneid_ti.rda")
load("rdata/cpm_2_geneid_cl.rda")
load("rdata/cpm_2_geneid_dcl.rda")

#################### PCA: log CPM > 2 matrix  #################### 
# Samples are in rows, transpose data frame
# Use Log2(CPM>2+1)
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
cpm_2_pca_cl = pca_func(cpm_2_geneid_cl)
cpm_2_pca_dcl = pca_func(cpm_2_geneid_dcl)

# Merge metadata with matrix by data ID
cpm_2_pca_plus_metadata = join(cpm_2_pca, metadata)
cpm_2_pca_ti_plus_metadata = join(cpm_2_pca_ti, metadata_ti)
cpm_2_pca_cl_plus_metadata = join(cpm_2_pca_cl, metadata_cl)
cpm_2_pca_dcl_plus_metadata = join(cpm_2_pca_dcl, metadata_dcl)

# Sanity check
table(cpm_2_pca_plus_metadata$SampleType)
table(cpm_2_pca_ti_plus_metadata$SampleType)
table(cpm_2_pca_cl_plus_metadata$SampleType)
table(cpm_2_pca_dcl_plus_metadata$SampleType)


#################### PCA: plot sample #################### 
pca_plot_sample <- function(cpm_2_pca_and_metadata) {
  # Get percent variance explained for plotting
  percentage <- round(cpm_2_pca_and_metadata$var_explained,2)
  percentage_formatted <- paste0("PC", paste(as.character(1:length(cpm_2_pca_and_metadata$var_explained))), " (", paste(as.character(percentage)), "%)")
  # Define the number of colors you want
  nb.cols <- length(unique(cpm_2_pca_and_metadata$SampleType))
  mycolors <- colorRampPalette(brewer.pal(3, "Set1"))(nb.cols)
  p_sample <- ggplot(cpm_2_pca_and_metadata,aes(x=PC1,y=PC2,color=SampleType)) + geom_point(size = 2) + theme_bw() + 
    xlab(percentage_formatted[1]) + ylab(percentage_formatted[2]) +
    ggtitle("PCA of 152 microRNA-seq libraries\nlabeled by sample type")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14, face="bold"),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  + scale_color_manual(values=mycolors)
  return(p_sample)
  
}

png(file = "figures/pca_log2cpmGreaterthan2_allSamples_SampleType.png",
    width = 10, 
    height = 8, units = "in", res = 1000)
pca_plot_sample(cpm_2_pca_plus_metadata)
dev.off()

#################### PCA: plot broader tissue type #################### 
pca_plot_broaderTissue <- function(cpm_2_pca_and_metadata) {
  # Get percent variance explained for plotting
  percentage <- round(cpm_2_pca_and_metadata$var_explained,2)
  percentage_formatted <- paste0("PC", paste(as.character(1:length(cpm_2_pca_and_metadata$var_explained))), " (", paste(as.character(percentage)), "%)")
  # Define the number of colors you want
  nb.cols <- length(unique(cpm_2_pca_and_metadata$BroaderTissueType))
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
  p_broadTissuetype <- ggplot(cpm_2_pca_and_metadata,aes(x=PC1,y=PC2,color=BroaderTissueType)) + geom_point(size = 3) + theme_bw() + 
    xlab(percentage_formatted[1]) + ylab(percentage_formatted[2]) +
   # Comment out when plotting all 152 
    # geom_label_repel(aes(label = BroaderTissueType),
    #                  box.padding   = 0.35,
    #                  point.padding = 0.5,
    #                  segment.color = 'grey50',
    #                  show.legend=FALSE) +
    
    ggtitle("PCA of 152 microRNA-seq libraries\nlabeled by derived tissue type")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  + scale_color_manual(values=mycolors)
  return(p_broadTissuetype)

}
png(file = "figures/pca_log2cpmGreaterthan2_TissueType_allSamples.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
pca_plot_broaderTissue(cpm_2_pca_plus_metadata)
dev.off()

png(file = "figures/pca_log2cpmGreaterthan2_TissueType_Tissues.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
pca_plot_broaderTissue(cpm_2_pca_ti_plus_metadata)
dev.off()

png(file = "figures/pca_log2cpmGreaterthan2_TissueType_CellLines.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
pca_plot_broaderTissue(cpm_2_pca_cl_plus_metadata)
dev.off()

png(file = "figures/pca_log2cpmGreaterthan2_TissueType_DifferentiatedCLs.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
pca_plot_broaderTissue(cpm_2_pca_dcl_plus_metadata)
dev.off()
########
########