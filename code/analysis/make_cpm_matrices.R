# 152 libraries across 102 experiments
# ENCODE 3 (34) and ENCODE 4 (68) experiments including tissue, differentiated cells, and cell lines
library(data.table)
setwd("/Users/liz/Documents/miRNA/ee282")

#################### Load counts files #################### 
load("rdata/counts_full_geneIDs.rda")
load("rdata/counts_full_geneids_tissues.rda")
load("rdata/counts_full_geneids_CLs_DCLs.rda")

#################### Convert counts matrix to CPM matrix #################### 
# Use nice function from google to convert to CPM (counts per million)
cpm_func <- function(counts) {
  cpm <- apply(counts,2, function(x) (x/sum(x))*1000000)
  cpm = as.matrix(cpm)
  rownames(cpm) = rownames(counts)
  colnames(cpm) = colnames(counts)
  return(cpm)
}

cpm_full_geneid = cpm_func(counts_geneid)
cpm_full_geneid_ti = cpm_func(counts_full_geneid_ti)
cpm_full_geneid_cl_dcl = cpm_func(counts_full_geneid_cl_dcl)

#################### Filter CPM matrices to CPM > 2 (at least one miRNA must be >2 CPM in a sample to be kept) #################### 
cpm_greaterthan2_func <- function(cpm) {
  check_max_cpm = max(cpm[1,])
  for (i in 2:dim(cpm)[1]) {
    check_max_cpm = c(check_max_cpm, max(cpm[i,]))
  }
  cpm_2 = cpm[which(check_max_cpm>2),]
  return(cpm_2)
}

cpm_2_geneid = cpm_greaterthan2_func(cpm_full_geneid) # 1055 miRNAs left
cpm_2_geneid_ti  = cpm_greaterthan2_func(cpm_full_geneid_ti) # 821 miRNAs left
cpm_2_geneid_cl_dcl = cpm_greaterthan2_func(cpm_full_geneid_cl_dcl) # 980 miRNAs left

save(cpm_2_geneid, metadata, file = "rdata/cpm_2_geneid.rda")
save(cpm_2_geneid_ti, metadata_ti, file = "rdata/cpm_2_geneid_ti.rda")
save(cpm_2_geneid_cl_dcl, metadata_cl_dcl, file = "rdata/cpm_2_geneid_cl_dcl.rda")

write.csv(cpm_2_geneid, file = "rdata/cpm_2_geneid.csv", quote = FALSE)
write.csv(cpm_2_geneid_ti, file = "rdata/cpm_2_geneid_ti.csv", quote = FALSE)
write.csv(cpm_2_geneid_cl_dcl, file = "rdata/cpm_2_geneid_cl_dcl.csv", quote = FALSE)
####
####