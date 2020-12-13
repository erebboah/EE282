library(Seurat)
setwd("/Users/liz/Documents/miRNA/ee282")
#################### UMAPs  #################### 
load("rdata/counts_0_geneIDs.rda")
table(colnames(counts_0_geneid) == metadata$miRNA_ENCODE_Data_ID) # sanity check

mirna <- CreateSeuratObject(counts = counts_0_geneid, project = "miRNA", min.cells = 0, min.features = 0)
mirna
mirna@meta.data$BroaderTissueType <- metadata$BroaderTissueType
mirna@meta.data$TissueType <- metadata$TissueType
mirna@meta.data$SampleType <- metadata$SampleType
mirna <- NormalizeData(mirna)
mirna <- FindVariableFeatures(mirna, selection.method = "vst", nfeatures = 100)
all.genes <- rownames(mirna)
mirna <- ScaleData(mirna, features = all.genes)
mirna <- RunPCA(mirna, features = VariableFeatures(object = mirna), npcs = 50)
ElbowPlot(mirna, ndims  = 50)
mirna <- RunUMAP(mirna, dims = 1:20)


nb.cols <- length(unique(mirna$BroaderTissueType))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)


png(file = "figures/umap_TissueType_allSamples.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
DimPlot(mirna, reduction = "umap", pt.size = 2, group.by = "BroaderTissueType", label=T, repel = T, cols = mycolors)+
            ggtitle("UMAP of 152 microRNA-seq libraries labeled by derived tissue type")+
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14),
                  plot.title = element_text(size = 14),
                  legend.title=element_text(size=14), 
                  legend.text=element_text(size=12),
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank())
dev.off()

png(file = "figures/umap_SampleType_allSamples.png",
    width = 10, 
    height = 8,units = "in", res = 1000)
DimPlot(mirna, reduction = "umap", pt.size = 2, 
        group.by = "SampleType",
        label=F, cols = c("#e41a19","#357eb7","#4cae48")) + 
    ggtitle("UMAP of 152 microRNA-seq libraries\nlabeled by sample type")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  
dev.off()


################ TISSUES only ################ 
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
mirna <- RunPCA(mirna, features = VariableFeatures(object = mirna), npcs = 50)
ElbowPlot(mirna, ndims  = 50)
mirna <- RunUMAP(mirna, dims = 1:20)

nb.cols <- length(unique(mirna$BroaderTissueType))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)


png(file = "figures/umap_TissueType_Tissues.png",
    width = 6, 
    height = 5,units = "in", res = 1000)
DimPlot(mirna, reduction = "umap", pt.size = 2, repel = T,group.by = "BroaderTissueType", label=T, cols = mycolors)+ 
    ggtitle("UMAP of 86 tissue microRNA-seq libraries\nlabeled by tissue type")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 14),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=12),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
dev.off()

#################### 
#################### 