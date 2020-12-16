#!/usr/bin/env Rscript
library("optparse",warn.conflicts = F)
library(ggplot2)

option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="Sequence statistics text file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=".", 
              help="Desired output folder [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="MicroRNA number", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

############## Plot ##############
nreads = as.data.frame(read.table(file=opt$data, sep = '\t',header=FALSE, stringsAsFactors = F))
colnames(nreads) = c("total","reads_with_bindingsite")
nreads$reads_with_bindingsite = gsub(" ","",nreads$reads_with_bindingsite)
nreads$nreads_with_bindingsite = as.numeric(sapply(strsplit((nreads$reads_with_bindingsite ), "/"), "[[", 1))
nreads$bamfile = sapply(strsplit((nreads$reads_with_bindingsite ), "intersect_"), "[[", 2)
nreads$bamfile = sapply(strsplit((nreads$bamfile ), ".bam.bed"), "[[", 1)

nreads$pct_reads = (nreads$nreads_with_bindingsite/nreads$total)*100

metadata = read.csv('/Users/liz/Documents/miRNA/ee282/ref/metadata_oct2020.csv', head=TRUE, sep=',', stringsAsFactors = FALSE)
metadata_filt = metadata[metadata$PB_ENCODE_BAM_ID %in% nreads$bamfile,]
metadata_filt2 = metadata_filt[!duplicated(metadata_filt$PB_ENCODE_BAM_ID), ]

pct_metadata  = merge(metadata_filt2, nreads, by.x = "PB_ENCODE_BAM_ID", by.y = "bamfile")

fname=paste0(opt$output,"/",opt$name,"_pct_read_barplot.png")

pct_metadata$PB_ENCODE_BAM_ID = factor(pct_metadata$PB_ENCODE_BAM_ID, levels=pct_metadata$PB_ENCODE_BAM_ID)
png(file = fname,
    width = 5, 
    height = 5, units = "in", res = 1000)
ggplot(pct_metadata, aes(x=PB_ENCODE_BAM_ID , y=pct_reads)) + ylab("% Reads") + 
  geom_bar(stat="identity", fill="#019f73") + theme(plot.margin=unit(c(0.2,0.1,0.1,0.8),"cm"))+
  ggtitle(paste0("Reads containing binding sites for miR-",opt$name,"\nin long read RNA-seq libraries\n(tissues only)"))+
  scale_x_discrete(labels=c(pct_metadata$TissueType)) +
  theme(legend.title=element_text(size=13),legend.text=element_text(size=12)) +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size=12)) + 
  theme(axis.text.y = element_text(size=12),axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1, size=12)) +
  theme(plot.title = element_text(size = 14, face="bold"))
dev.off()

#####
