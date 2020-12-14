# 86 tissue samples
# DE of heart tissue compared to all other tissue samples

library(ggplot2)
library(EnhancedVolcano)
library(edgeR)

setwd("/Users/liz/Documents/miRNA/ee282")

#################### EdgeR - simple exact test between heart tissue and all others #################### 
# EdgeR takes in counts
load("rdata/counts_0_geneids_tissues.rda")

metadata_ti$Heart[metadata_ti$BroaderTissueType == "Heart"]<- "Heart"
metadata_ti$Heart[metadata_ti$BroaderTissueType != "Heart"]<- "NotHeart"

group <- metadata_ti$Heart
y <- DGEList(counts=counts_0_geneid_ti, group=group)
y$samples

# Filter
keep <- filterByExpr(y, group=group)
table(keep) # keep 637 miRNAs
y <- y[keep, , keep.lib.sizes=FALSE]

# Internal edgeR normalization
y <- calcNormFactors(y)

# Classic edgeR pipeline: pairwise comparison between 2 groups
y <- estimateDisp(y)

et <- exactTest(y, pair = c("NotHeart","Heart"))
topTags(et)
edger_result_heart = et$table

#################### EdgeR - simple exact test between stem cells and all other cell lines #################### 
# EdgeR takes in counts
load("rdata/counts_0_geneids_cls_dcls.rda")

metadata_cl_dcl$SC[metadata_cl_dcl$BroaderTissueType == "hESC" | metadata_cl_dcl$BroaderTissueType == "iPSC"] <- "SC"
metadata_cl_dcl$SC[metadata_cl_dcl$BroaderTissueType != "SC"]<- "NotSC"

group <- metadata_cl_dcl$SC
y <- DGEList(counts=counts_0_geneid_cl_dcl, group=group)
y$samples

# Filter
keep <- filterByExpr(y, group=group)
table(keep) # keep 587 miRNAs
y <- y[keep, , keep.lib.sizes=FALSE]

# Internal edgeR normalization
y <- calcNormFactors(y)

# Classic edgeR pipeline: pairwise comparison between 2 groups
y <- estimateDisp(y)

et <- exactTest(y, pair = c("NotSC","SC"))
topTags(et)
edger_result_sc = et$table
  
  
#################### EdgeR - heart tissue vs all tissues; stem cells vs all lines  #################### 
p1 = EnhancedVolcano(edger_result_heart,
                    lab = rownames(edger_result_heart),
                    title = 'microRNAs upregulated (positive LFC)\nin heart tissue',
                    subtitle = '-Log10(p-val) < 10, |Log2FC| > 2\n354 microRNAs total',
                    x = 'logFC',
                    y = 'PValue',
                    FCcutoff = 2,
                    pCutoff = 10e-11,
                    boxedLabels = TRUE,
                col = c("#717171", "#4daf4a", "#377eb7", "#e41a1b"),
                pointSize = 3.0,
                labSize = 5.0,
                legendLabSize = 13,
                titleLabSize = 17,
                subtitleLabSize = 13,
                legendPosition = "bottom",
                caption = NULL,
                captionLabSize = 13,
                    drawConnectors = TRUE,
                    colAlpha = 1)

p2 = EnhancedVolcano(edger_result_sc,
                     lab = rownames(edger_result_sc),
                     title = 'microRNAs upregulated (positive LFC)\nin stem cells',
                     subtitle = '-Log10(p-val) < 12, |Log2FC| > 2\n587 microRNAs total',
                     x = 'logFC',
                     y = 'PValue',
                     FCcutoff = 2,
                     pCutoff = 10e-13,
                     boxedLabels = TRUE,
                     col = c("#717171", "#4daf4a", "#377eb7", "#e41a1b"),
                     pointSize = 3.0,
                     labSize = 5.0,
                     legendLabSize = 13,
                     titleLabSize = 17,
                     subtitleLabSize = 13,
                     legendPosition = "bottom",
                     caption = NULL,
                     captionLabSize = 13,
                     drawConnectors = TRUE,
                     colAlpha = 1)

png(file = "figures/volcano_edgeR_heart_stemCells.png",
    width = 12, 
    height = 8, units = "in", res = 1000)
ggarrange(p1,p2,common.legend = TRUE,legend="bottom")
dev.off()

#################### Save results #################### 
save(edger_result_heart, file="rdata/edger_result_heartTI_vs_allotherTI.csv",quote=F)
save(edger_result_sc, file="rdata/edger_result_stemCells_vs_allotherCLDCL.rda", quote=F)


#####
#####
#####
