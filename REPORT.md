##  Final Project Report
#### Elisabeth Rebboah
#### Eco Evo 282
#### December 16, 2020
***
## Methods
### MicroRNA
MicroRNA-seq is a form of RNA-seq that captures small, functional RNAs using a 3' adaptor that binds specifically to their 3' hydroxyl group, as oppposed to mRNA-seq that often relies on a 3' poly-A tail.

![Figure 1](fig1_experimentOverview.png)

The libraries are processed by trimming the long adapters on the 3' and 5' ends using cutadapt, leaving around 22 base pairs to be mapped to the genome and quantified using STAR. The ENCODE quality control standards require at least 5 million mapped/multi-mapped reads, at least 300 microRNAs expressed over 2 counts per million (CPM), and a Spearman correlation > 0.85 if there are biological replicates.

The tsv files containing quantifications for 152 human microRNA-seq libraries available on the [ENCODE website](https://www.encodeproject.org/matrix/?type=Experiment&status=released&perturbed=false&assay_title=microRNA-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.rfa=ENCODE3&award.rfa=ENCODE4&perturbed=true&status=submitted) were downloaded. A detailed metadata file used throughout this project found in the data folder was used to combine the libraries into one counts matrix of 152 samples by 1,881 microRNAs (hg38 GENCODE v29). Counts were also converted to CPM and filtered such that each microRNA had at least 2 CPM in at least one sample, leaving 1,055 microRNAs. The scripts used to generate counts matrices and CPM matrices can be found in code/analysis and are called make_counts_matrices.R and make_cpm_matrices.R, respectively, and the R data/csv files can be found in data/processed. 

Principal component analyis (PCA) was run on this CPM>2 matrix containing all samples, as well as on tissues only (86 samples) and cell lines/differentiated cell lines only (66 samples). In order to represent more variance and have a better understanding of the clusters, a non-linear dimensionality reduction method called Uniform Manifold Approximation and Projection (UMAP) was also implemented on the counts matrices via the Seurat workflow. This workflow normalizes the data, hence why raw counts were used. The script used to generate PCA and UMAP plots can be found in code/analysis and is called pca_umap_figures.R.     

Based on PCA and UMAP analysis, differential expression (DE) analysis using edgeR was carried out on heart tissue samples compared to all other tissues as well as stem cells compared to all other cell line/differentiated cell line samples. Different p-value/log fold change cutoffs were used in the two volcano plots to be visually appealing; full edgeR results can be found in data/processed. 

### Long read RNA
First, the long read RNA-seq libraries that match (were built from the same RNA) as microRNA-seq libraries were found using the metadata file. The building/submission of long-read libraries and microRNA-seq libraries is an ongoing process, so there are many samples where the long-read libraries exist or are on schedule to be sequenced and they are labeled as "waiting" in the simple stacked bar plots made from the metadata file using the script code/analysis/stacked_barplot_matchingPB.R. Bam files for the 23 matching samples (11 tissues, 12 cell lines/diff. cell lines) were downloaded from the ENCODE website.

Next, the bash script code/scripts/find_binding_sites_from_bam_TISSUES.sh was used to find long reads containing a microRNA binding site for a microRNA of interest. The script uses awk to filter the predicted targets database file (data/Predicted_Target_Locations.default_predictions.hg38.named.bed) for a specific microRNA, then loops through the matching tissue bam files, checking the binding site bed file against each one using bedtools intersect. The script outputs the filtered database file, the number of reads per library containing a binding site and total number of reads, and bed files containing the genomic location of the binding site, the long read name, and the genomic location of the long read. The bash script calls an R script to a simple bar plot of the fraction of reads with a specific microRNA binding site per long read RNA-seq library. The process is the same for code/scripts/find_binding_sites_from_bam_LINES.sh, except the matching libraries from cell line/differentiated cell lines were used instead of matching tissue libraries. 

## Results
PCA of all 152 samples shows that samples mainly cluster together based on if the RNA was isolated from from tissues versus cell lines and differentiated cell lines. This could be caused by differing RNA quality between samples; due to the experimental processing, microRNA-seq performs better on high-quality (RIN>8) RNA. In comparison, the UMAP clusters are more tightly grouped by tissue type. For example, the UMAP shows the myocyte differentiated cells cluster closely to muscle and heart tissues, while this distinction is less clear in the PCA.

![Figure 2 PCA and UMAP of all samples](pca_umap_allSamples_SampleTypes.png)
![Figure 3 PCA and UMAP of all samples labeled by tissue](pca_umap_allSamples_TissueTypes.png)
![Figure 4 PCA and UMAP of tissue samples](pca_umap_TissueType_Tissue.png)
![Figure 5 PCA and UMAP of cell line/differentiated cell line samples](pca_umap_TissueType_CLsDCLs.png)

The differentially expressed microRNAs are indicated in the following volcano plots. Using the same p-value cutoff of 0.01 and a log fold change cutoff of 2, there were 57 differentially expressed microRNAs in the heart vs. all tisses and 197 differentially expressed microRNAs in stem cells vs. all cell lines/diff. cell lines. 

![Figure 6 DE analysis with edgeR](volcano_edgeR_heart_stemCells.png)

There are 11 matching tissue long-read libraries, mostly from heart, and 12 matching cell line long read libraries, somewhat scattered but there are 4 matching stem cell libraries. To note, the metadata file was last curated in October 2020 and there have been both microRNA and long read library submissions since then.

![Figure 7 Matching long read libraries](matchedPB_barplots.png)

## Discussion
As an initial exploratory analysis, all microRNA-seq samples were projected in low-dimensionality spaces. MicroRNAs expressed in samples from similar derived/simulated tissues that clustered together and apart from other samples are likely to be tissue-specific. The central role of muscle-specific microRNAs could be why the muscle samples form their own cluster. In addition, microRNAs have been implicated to have a role in stem cell self-renewal, possibly reflected by the hESC/iPSC samples clustering together. 

DE analysis revealed several known microRNAs for the tissue of interest, both upregulated and downregulated. For example, miR-208A  is encoded within an intron of a muscle myosin heavy chain gene, Myh6. Found as one of the top downregulated microRNAs in heart, miR-200C overexpression inhibits myogenic differentiation. In stem cells, the upregulated microRNAs miR-518b and miR-520g have support in literature, and accumulation of the downregulated microRNA let-7 is prevented by Lin-28 which promotes pluripotency.

From the DE results and literature search, a few microRNAs were chosen to investigate further.
