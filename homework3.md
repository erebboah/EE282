# Elisabeth Rebboah
# Eco Evo 282
# Homework 3
***
## Summarize Genome Assembly
### File Integrity
The command md5sum can be used to check the integrity of downloaded files by comparing it to the md5sums listed on the FlyBase website. Saving the md5sum of the file as a variable makes it easier to search for the exact md5sum sequence in the reference list of md5sums from the website without making a mistake copy-pasting.
```
md5=($(md5sum dmel-all-chromosome-r6.36.fasta.gz))
grep $md5 md5sum.txt
```
### Calculate Summaries of the Genome
```
faSize dmel-all-chromosome-r6.36.fasta.gz
Total number of nucleotides, Ns, and sequences:
143726002 bases (1152978 Ns 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real
```
1. From the convenient output of faSize, there are 143,726,002 nucleotides in the genome.  
2. There are 1,152,978 Ns.
3. There are 1,870 total sequences.

The files can be downloaded and above code can be run in 1 line, outputting a folder called fly_ref in the same directory:
```
bash hw3_genome_summary.sh
```
## Summarize an Annotation File
### File Integrity
As with the genome, the annotation file integrity can be verified with the same code:
```
md5=($(md5sum dmel-all-r6.36.gtf.gz))
grep $md5 md5sum.txt
```
The previous script renames the genome md5sum.txt to md5sum_genome.txt to avoid confusion.

### Compile a Report Summarizing the Annotation
1. The GTF must first be unzipped, then the file can be cut to only show the third field, which contains features. Next the features are sorted, which outputs a long list of ordered features, and uniq -c counts the number of occurences of each feature. Finally, sort -nr sorts numerically and in reverse in order to display results from most to least common.
```
cut -f3 dmel-all-r6.36.gtf | sort | uniq -c | sort -nr
 189268 exon
 162578 CDS
  46664 5UTR
  33629 3UTR
  30812 start_codon
  30754 stop_codon
  30728 mRNA
  17875 gene
   3047 ncRNA
    485 miRNA
    366 pseudogene
    312 tRNA
    300 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
```
2. The chromosome column and the feature columns are pulled out then filtered to only the gene features using grep, making sure that the feature starts with the word gene, unlike pseudogene. Next, grep is used again for each of the 7 chromosomes, printing the matched lines. Finally, as above, sort followed by uniq -c counts the number of occurences for each chromosome. 
```
cut -f1,3 dmel-all-r6.36.gtf | grep -E '\<gene' | grep -Eo 'X|Y|[2,3][L,R]|4' | sort | uniq -c
   3516 2L
   3653 2R
   3486 3L
   4225 3R
    125 4
   2691 X
    113 Y
```

The files can be downloaded and above code can be run in 1 line, adding to the folder called fly_ref in the same directory:
```
bash hw3_annotation_summary.sh
```
