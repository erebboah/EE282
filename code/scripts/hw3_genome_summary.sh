#!/bin/bash

source activate ee282
mkdir fly_ref/
cd fly_ref/

### Genome File Integrity
echo Downloading md5sum information and gzipped fasta for all chromosomes...
# Get md5sum information from FlyBase
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/md5sum.txt

# Download gzipped fasta for all chromosomes
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz

echo Make sure checksum of downloaded file matches reference:
md5=($(md5sum dmel-all-chromosome-r6.36.fasta.gz))
grep $md5 md5sum.txt
mv md5sum.txt md5sum_genome.txt

### Genome Summary
echo Total number of nucleotides, Ns, and sequences:
faSize dmel-all-chromosome-r6.36.fasta.gz
