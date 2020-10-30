#### Elisabeth Rebboah
#### Eco Evo 282
#### Project Analysis Proposal
***
##### Datasets
I will be using the database of 86 microRNA-seq libraries built from human tissue samples available on the [ENCODE website.](https://www.encodeproject.org/matrix/?type=Experiment&status=released&perturbed=false&assay_title=microRNA-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.rfa=ENCODE3&award.rfa=ENCODE4). Although interesting questions could be asked of the differences between microRNA expression in cell line and tissue samples, a more focused analysis on tissue samples only could reveal known and novel tissue-specific microRNAs rather than technical differences between _in vivo_ and _in vitro_ samples.

In addition, I will be using 11 long read RNA-seq libraries built from some of the same tissues as the microRNA-seq samples, including: 2 heart ventricle samples, one vena cava, one aorta, 3 lung samples, colon mucosa, mesenteric fat pad, adrenal gland, and ovary. A metadata file contains biosample information, data IDs, experiment IDs, etc. for both microRNA and long read data.

##### Analysis/Figures
First the quantification tsv files from each library (for example, [ENCFF428XME](https://www.encodeproject.org/files/ENCFF428XME/)) will be concatenated into a single matrix using R. The counts matrix for the 86 tissue samples will be analyzed using the R edgeR differential expression package to compare heart samples to all other tissues. The resulting differentially expressed microRNAs will be plotted on a volcano plot using the R package EnhancedVolcano. 

From this list of microRNAs upregulated in cardiac tisue, I will choose several that have binding sites in the ["Default predictions" BED file.](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) The BED file has to be converted from hg19 to hg38 using [LiftOver.](https://genome.ucsc.edu/cgi-bin/hgLiftOver)

In order to count the number of transcripts in each long read library that include at least one binding site, I will first use grep to subset the hg38 BED file for a microRNA of interest: 
```console
#!/bin/bash

mirna=$1

grep ${mirna} /Users/liz/Documents/miRNA/ee282/ref/Predicted_Target_Locations.default_predictions.hg38.named.bed > /Users/liz/Documents/miRNA/ee282/ref/Predicted_Target_Locations.default_predictions.hg38.named.${mirna}.bed 
```
Next I will use bedtools intersect to subset each GTF for the matching long reads (for example, [ENCFF010KII](https://www.encodeproject.org/files/ENCFF010KII/)).
```console
cd /Users/liz/Documents/miRNA/ee282/gtf_files/tissues/
for i in `ls -1v *.gtf`; do
        bedtools intersect -wa -F 1 -a ${i%} -b /Users/liz/Documents/miRNA/ee282/ref/Predicted_Target_Locations.default_predictions.hg38.named.${mirna}.bed > /Users/liz/Documents/miRNA/ee282/output/${i%} 
done
```
And finally I can subset the GTF with awk to get columns of interest, such as TALON transcript ID which I can cross-reference with transcript quantifications, such as  [ENCFF304FFP](https://www.encodeproject.org/files/ENCFF304FFP/). 
```console
awk -F "\t" '$3 == "transcript" { print $9 }' ${i%}  > ${i%}_transcripts.txt
awk -F "talon_transcript" '{ print $2 }' ${i%}_transcripts.txt | awk -F "; " '{print $1}' > talontranscriptID
```
Using wc -l, I will count the number of transcripts in each of the 11 samples that contain a binding site compared to the total number of transcripts detected in the sample. I will make bar plots in R using ggplot showing the fraction of "known" or "novel" transcripts of the ones containing a binding site for each of the 11 libraries and include the percentage of transcripts with a binding site out of the total number detected.

##### Conclusion
I expect heart samples to be enriched in binding sites for microRNAs upregulated in cardiac tissue, but the results may be unclear because of the quality of the predicted targets database. Regardless, this preliminary analysis is a good start at integrating microRNA-seq and long reads with the potential to be built upon in the future using a different database and/or a different format of long read quantifications rather than GTF.
