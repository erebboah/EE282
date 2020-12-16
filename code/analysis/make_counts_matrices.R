# 152 libraries across 102 experiments
# ENCODE 3 (34) and ENCODE 4 (68) experiments including tissue, differentiated cells, and cell lines
library(data.table)
setwd("/Users/liz/Documents/miRNA/ee282")

#################### Load metadata #################### 
metadata = read.csv('ref/metadata_oct2020.csv', head=TRUE, sep=',', stringsAsFactors = FALSE)

#################### Read in .tsv and make counts matrix across all samples in data folder #################### 
samples=metadata$miRNA_ENCODE_Data_ID

thisMat=read.table(paste0("/Users/liz/Documents/miRNA/ENCODE/data/",samples[1],".tsv"),sep="",head=FALSE)
counts=thisMat[-c(1:4),c(1,2)]

for(i in 2:length(samples)){
  thisMat=read.table(paste0("/Users/liz/Documents/miRNA/ENCODE/data/",samples[i],".tsv"),sep="",head=FALSE)
  thisMat=thisMat[na.omit(match(counts$V1,thisMat$V1)),]
  counts=as.data.frame(cbind(counts,thisMat[,"V2"]))
}

dim(counts)
# Set rownames as ENSEMBL IDs
rownames(counts) = counts$V1
counts$V1 = NULL
# Set rownames as data IDs
colnames(counts) = metadata$miRNA_ENCODE_Data_ID
head(counts) 
dim(counts) # 1,881 miRNA by 152 samples

# Save full counts matrix with ensembl IDs
save(counts, metadata, file="rdata/counts_full_ensemblIDs.rda")
write.csv(counts, file="rdata/counts_full_ensemblIDs.csv", quote = F)

# Use metadata to subset counts into tissues, cell lines, and differentiated cell lines
metadata_ti = metadata[metadata$SampleType == "Tissue",]
counts_full_ensembl_ti = counts[,colnames(counts) %in% metadata_ti$miRNA_ENCODE_Data_ID]
metadata_ti_sorted=metadata_ti[order(match(metadata_ti$miRNA_ENCODE_Data_ID,colnames(counts_full_ensembl_ti))),]
table(metadata_ti$miRNA_ENCODE_Data_ID == colnames(counts_full_ensembl_ti))
save(counts_full_ensembl_ti, metadata_ti, file="rdata/counts_full_ensemblIDs_tissues.rda")
write.csv(counts_full_ensembl_ti, file="rdata/counts_full_ensemblIDs_tissues.csv", quote = F)

metadata_cl_dcl = metadata[metadata$SampleType == "CellLine" | metadata$SampleType == "DifferentiatedCells",]
counts_full_ensembl_cl_dcl = counts[,colnames(counts) %in% metadata_cl_dcl$miRNA_ENCODE_Data_ID]
metadata_ti_sorted=metadata_cl_dcl[order(match(metadata_cl_dcl$miRNA_ENCODE_Data_ID,colnames(counts_full_ensembl_cl_dcl))),]
table(metadata_cl_dcl$miRNA_ENCODE_Data_ID == colnames(counts_full_ensembl_cl_dcl))
save(counts_full_ensembl_cl_dcl, metadata_cl_dcl, file="rdata/counts_full_ensemblIDs_CLs_DCLs.rda")
write.csv(counts_full_ensembl_cl_dcl, file="rdata/counts_full_ensemblIDs_CLs_DCLs.csv", quote = F)


#################### Change ensembl ID to gene names in counts matrices ####################
# read in annotation file make from mapping GTF 
annot  = as.data.frame(read.table('ref/miRNA_v29_ENSG_genename.gtf', head=F))
annot$V2=NULL # remove semicolon
annot$V4=NULL # remove semicolon
# annot is ENSG and matching mature miRNA ID
colnames(annot) = c("ensembl","geneid")
dim(annot)
head(annot)

# Change IDs in counts matrix
# Make a column of ensembl IDs
counts$ensembl = rownames(counts)
# merge annotation with counts matrix to get matching miRNA ID
counts_geneid <- merge(counts, annot, by='ensembl')
geneid_rownames = counts_geneid$geneid
# Get rid of extra columns 
counts_geneid$ensembl <- NULL
counts_geneid$geneid <- NULL
# Have to make it a matrix to have duplicate rownames
counts_geneid = as.matrix(counts_geneid, quote=F)
rownames(counts_geneid) = geneid_rownames
dim(counts_geneid)
head(counts_geneid)

save(counts_geneid, metadata, file="rdata/counts_full_geneIDs.rda")
write.csv(counts_geneid, file="rdata/counts_full_geneIDs.csv", quote = F)

counts_full_geneid_ti = counts_geneid[,colnames(counts_geneid) %in% metadata_ti$miRNA_ENCODE_Data_ID]
counts_full_geneid_cl_dcl = counts_geneid[,colnames(counts_geneid) %in% metadata_cl_dcl$miRNA_ENCODE_Data_ID]

# Save gene name annotated counts matrices
save(counts_full_geneid_ti,metadata_ti, file="rdata/counts_full_geneids_tissues.rda")
save(counts_full_geneid_cl_dcl, metadata_cl_dcl, file="rdata/counts_full_geneids_CLs_DCLs.rda")

write.csv(counts_full_geneid_ti, file="rdata/counts_full_geneids_tissues.csv", quote = F)
write.csv(counts_full_geneid_cl_dcl, file="rdata/counts_full_geneids_CLs_DCLs.csv", quote = F)

#################### Filter counts matrixes matrices to counts > 0 #################### 
counts_greaterthan0_func <- function(counts) {
  counts_0 = counts[rowSums(counts[,-1])>0,]
  return(counts_0)
}

counts_0_geneid = counts_greaterthan0_func(counts_geneid)
counts_0_geneid_ti = counts_greaterthan0_func(counts_full_geneid_ti)
counts_0_geneid_cl_dcl = counts_greaterthan0_func(counts_full_geneid_cl_dcl)

# Save gene name annotated counts matrices
save(counts_0_geneid, metadata, file="rdata/counts_0_geneIDs.rda")
save(counts_0_geneid_ti,metadata_ti, file="rdata/counts_0_geneids_tissues.rda")
save(counts_0_geneid_cl_dcl, metadata_cl_dcl, file="rdata/counts_0_geneids_CLs_DCLs.rda")

write.csv(counts_0_geneid, file="rdata/counts_0_geneIDs.csv", quote = F)
write.csv(counts_0_geneid_ti, file="rdata/counts_0_geneids_tissues.csv", quote = F)
write.csv(counts_0_geneid_cl_dcl, file="rdata/counts_0_geneids_CLs_DCLs.csv", quote = F)

#################### 
#################### 
#################### 