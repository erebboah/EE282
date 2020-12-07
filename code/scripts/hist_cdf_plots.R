#!/usr/bin/env Rscript
library("optparse",warn.conflicts = F)
library(ggplot2)
library(ggpubr)

option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="Sequence statistics text file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=".", 
              help="Desired output folder [default= %default]", metavar="character"),
  make_option(c("-b", "--bins"), type="character", default=NULL, 
              help="Number of bins in histograms", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="Desired output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

############## Plot ##############
seq_stats = as.data.frame(read.table(file=opt$data, sep = '\t',header=FALSE, stringsAsFactors = F))
colnames(seq_stats) = c("Name","Length","GC")

p1 = ggplot(seq_stats, aes(x=log2(Length))) + 
  geom_histogram(color="black", fill="#57b4e9", bins = opt$bins) + 
  labs(title= paste0("Histogram of sequence length distribution\nfor sequences ", opt$name), x = "Log2(Sequence length)", y = "Number of sequences") +
  theme(text = element_text(size=14))

p2 = ggplot(seq_stats, aes(x=GC)) + 
  geom_histogram(color="black", fill="#e79f00", bins = opt$bins) + 
  labs(title= paste0("Histogram of sequence GC% distribution\nfor sequences ", opt$name), x = "GC%", y = "Number of sequences") +
  theme(text = element_text(size=14))
  
plots = ggarrange(p1, p2, ncol = 2, nrow = 1)

fname=paste0(opt$output,"/",opt$name,"_histograms.png")
ggsave(fname, plots, width = 10, height = 5)


#####


