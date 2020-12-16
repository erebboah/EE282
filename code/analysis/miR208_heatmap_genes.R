library(ggplot2)
library(stringr)
library(gplots)
library(RColorBrewer)
library(pheatmap)

setwd("/Users/liz/Documents/miRNA/ee282/")
metadata = read.csv('ref/metadata_oct2020.csv', head=TRUE, sep=',', stringsAsFactors = FALSE)

samples = dir("/Users/liz/Documents/miRNA/ee282/great/tissues")

great=read.table(paste0("great/tissues/",samples[1]))
counts = as.data.frame(table(great$V2))
name = sapply(strsplit(samples[1], "intersect_"), "[[", 2)
name = gsub("_great.txt","",name)
colnames(counts) = c("gene",name)

for(i in 2:length(samples)){
  great=read.table(paste0("great/tissues/",samples[i]))
  counts2 = as.data.frame(table(great$V2))
  name = sapply(strsplit(samples[i], "intersect_"), "[[", 2)
  name = gsub("_great.txt","",name)
  colnames(counts2) = c("gene",name)
  
  counts = merge(counts,counts2,by="gene", all=T)
}
rownames(counts) = counts$gene
counts$gene = NULL

dim(counts)
counts = na.omit(counts)

colnames(counts) = metadata[match(colnames(counts),metadata$PB_ENCODE_BAM_ID),"TissueType"]

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",
#                   ns=asNamespace("pheatmap"))

png(file = "figures/miR208_bindingsite_heatmap.png",
    width = 9, 
    height = 4, units = "in", res = 1000)
pheatmap(t(counts), scale="column",fontsize_row = 9, fontsize_col = 9,
         main = "Number of miR-208 binding sites per gene (scaled)",treeheight_col=25,treeheight_row=25)
dev.off()
###
###
