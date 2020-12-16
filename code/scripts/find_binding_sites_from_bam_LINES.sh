#!/bin/bash

mirna_number=$1

path="/Users/liz/Documents/miRNA/ee282/"

# Search 4th column of bed file (names field) for miRNA of interest
# Have to consider instances where multiple miRNAs have same target
# Output filtered bed file 
search1=/"$mirna_number"
search2=-"$mirna_number"
awk -v s1=$search1 -v s2=$search2 \
'$4 ~ s1 || $4 ~ s2' "$path"/ref/Predicted_Target_Locations.default_predictions.hg38.named.bed > "$path"/ref/Targets_miR-"$mirna_number".bed

# Bedtools intersect with specific miRNA binding sites and long read RNA-seq bam file
# Output bed file
cd "$path"/bam/cellLines
for i in `ls -1v *.bam`
	do bedtools intersect -bed -abam "$path"/bam/cellLines/${i%} \
	-b "$path"/ref/Targets_miR-"$mirna_number".bed | awk -F "\t" '{print $1 "\t" $7 "\t" $8}' > "$path"/output/cellLines/miR-"$mirna_number"_intersect_${i%}.bed
done

# Count number of reads in bam file and number of reads with binding site
# Output text file 
cd "$path"/bam/cellLines
num_reads_total=$(for i in `ls -1v *.bam`; do
          samtools view -c "$path"/bam/cellLines/${i%}
      done)

cd "$path"/bam/cellLines
num_reads_intersect=$(for i in `ls -1v *.bam`; do
          wc -l "$path"/output/cellLines/miR-"$mirna_number"_intersect_${i%}.bed
      done)

paste <(printf %s "$num_reads_total") <(printf %s "$num_reads_intersect") > "$path"/output/cellLines/miR-"$mirna_number"_NumReadOverlaps.txt

Rscript "$path"/scripts/percent_bindingSite_barplot_LINES.R -d "$path"/output/cellLines/miR-"$mirna_number"_NumReadOverlaps.txt -o "$path"/output/cellLines -n "$mirna_number"
