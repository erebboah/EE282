##  Final Project Report
#### Elisabeth Rebboah
#### Eco Evo 282
#### December 16, 2020
***
## Methods
MicroRNA-seq is a form of RNA-seq that captures small, functional RNAs using a 3' adaptor that binds specifically to their 3' hydroxyl group, as oppposed to mRNA-seq that often relies on a 3' poly-A tail.

![Figure 1](fig1_experimentOverview.png)

The libraries are processed by trimming the long adapters on the 3' and 5' ends using cutadapt, leaving around 22 base pairs to be mapped to the genome and quantified using STAR. The ENCODE quality control standards require at least 5 million mapped/multi-mapped reads, at least 300 miRNAs expressed over 2 counts per million (CPM), and a Spearman correlation > 0.85 if there are biological replicates.

The tsv files containing quantifications for all 152 human microRNA-seq libraries available on the [ENCODE website](https://www.encodeproject.org/matrix/?type=Experiment&status=released&perturbed=false&assay_title=microRNA-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.rfa=ENCODE3&award.rfa=ENCODE4&perturbed=true&status=submitted) were downloaded. A detailed metadata file used throughout this project was used to combine the libraries into one counts matrix of 152 samples by 1,881 microRNAs (hg38 GENCODE v29). Counts were also converted to CPM and filtered such that each microRNA had at least 2 CPM in at least one sample, leaving 1,055 microRNAs. The scripts used to generate counts matrices and CPM matrices can be found in code/scripts and are called make_counts_matrices.R and make_cpm_matrices.R, respectively. 

Principal component analyis (PCA) was run on this CPM>2 matrix, as well as on CPM>2 matrices for each sample type: cell line (36 samples), differentiatd cells (30 samples), and tissues (86 samples). In order to represent more variance and have a better understanding of the clusters, a non-linear dimensionality reduction method called Uniform Manifold Approximation and Projection (UMAP) was also implemented on the counts matrices via the Seurat workflow. The UMAP method normalizes the data, hence why counts were used. The scripts used to generate PCA and UMAP plots can be found in code/scripts and are called pca_figures.R and umap_figures.R, respectively.      


## Results
PCA of all 152 samples shows that samples mainly cluster together based on if the RNA was isolated from from tissues versus cell lines and differentiated cell lines. This could be due to differing RNA quality between samples; due to the experimental processing, microRNA-seq performs better on high-quality (RIN>8) RNA. In comparison, the UMAP clusters are more tightly grouped by tissue type. For example, the myocyte differentiated cells samples cluster closely to muscle and heart tissue samples, and pancreas tissue samples cluster closely to pancreatic differentiated cells.

## Discussion
As an initial exploratory analysis, all microRNA-seq samples were projected in low-dimensionality spaces. MicroRNAs expressed in samples from similar derived/simulated tissues that clustered together  and apart from other samples are likely to be tissue-specific.  entral role of muscle-specific microRNAs or myomiRs the muscle samples were as well as hormonal tissues such as pancreas, ovaries, testes, and adrenal glands, and stem cells. Based on PCA and UMAP analysis, differential expression analysis was carried out on cardiac tissue samples compared to all other tissues.  and are interesting candidates for further analysis.
