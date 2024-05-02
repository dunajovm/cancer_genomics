# Variant calling workflow

In this notebook workflow used for variant calling will be described. Pipeline will run mainly in batch 
command-line.

### Getting the files and the programs
Pipeline will be runing in its own virtual enviroment - to be sure, that none of the dependencies 
for this pipeline can create problems with already installed programs.
```bash
virtualenv exercise_1
source exercise_1/bin/activate
sudo apt install -y bwa
sudo apt install -y bcftools
```
Now we will get files, which we will be working with, locally. Namelly:
- reference file for chromozom 7 in fasta format - *ch7.fa.gz*
- forward reads on chromzom 7 in fastq format - *R1.fastq.gz*
- reverse reads on chromzom 7 in fastq format - *R2.fastq.gz*

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
wget https://gear-genomics.embl.de/data/.slides/R1.fastq.gz
wget https://gear-genomics.embl.de/data/.slides/R2.fastq.gz
```

### Aligning reads to reference using bwa
We will assume, that the reads are in good quality, in no need of trimming and filtering.
Therefore, now we will do mapping reads to the reference genome using bwa aligner.
Namely, we will use bwa-mem, as it is quick and accurate.

```bash

```
