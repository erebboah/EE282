#!/bin/bash

cd fly_ref/

### Annotation File Integrity
echo Downloading md5sum information and current GTF for D. melanogaster...
# Get md5sum information from FlyBase
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/md5sum.txt

# Download gzipped GTF
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/dmel-all-r6.36.gtf.gz

echo Make sure checksum of downloaded file matches reference:
md5=($(md5sum dmel-all-r6.36.gtf.gz))
grep $md5 md5sum.txt
mv md5sum.txt md5sum_annot.txt

### Annotation Summary
echo Total number of features of each type:
gunzip dmel-all-r6.36.gtf.gz
cut -f3 dmel-all-r6.36.gtf | sort | uniq -c | sort -nr

echo Total number of genes per chromosome arm:
cut -f1,3 dmel-all-r6.36.gtf | grep -E '\<gene' | grep -Eo 'X|Y|[2,3][L,R]|4' | sort | uniq -c
gzip dmel-all-r6.36.gtf
