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
sudo apt install -y samtools
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
Bwa-mem needs not only reference genome, but also the index file of the reference genome.
Index file will be generated via bwa-index.

```bash
bwa index chr7.fa.gz
bwa mem chr7.fa.gz R1.fastq.gz R2.fastq.gz > alignment.bam
```
### Variant calling using bcftools
We will do variant calling using bcftools variant caller.
Bcftools variant caller needs sorted alignment, therefore we will first sort the alignment using samtools sort.
For variant calling, we will asume, that the organism is diploid. We will generate variants in vcf format.
We will give caller the -v option, which will give us only variant sites in the output file.

```bash
samtools sort -o alignment.sorted.bam alignment.bam 
bcftools mpileup -f chr7.fa alignment.sorted.bam | bcftools call -vm -Ov -o variants.vcf

```
We can take a look at the resulting variants via:
```bash
bcftools view variants.vcf
# or
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" variants.vcf
# number of found variants is 41
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" variants.vcf | wc -l

```
### Anotating variants via VEP
We will use online version of VEP for anotating variants.
VEP can be found in following [link](https://www.ensembl.org/Tools/VEP).
As the reference, we will use human genome, assembly GRCh38.p14.


![image1](./images/image1.png)


Result showed, that one gene was affected by the found variants.


![image2](./images/image2.png)

Affected was gene [CARD11](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000198286;r=7:2906142-3044228;tl=B5j5OQNfJqj0ZkRK-10087084). Out of total 41 variants, 38 variants were located in non-coding intron or exon parts of the gene. Only 3 variants were located in the coding part of the protein transcripts. Out of them, 2 were synonimous variants, 1 was stop-gained variant.

![image3](./images/image3.png)

The result of the VEP annotation can be found in VCF format [here](./images/B5j5OQNfJqj0ZkRK.vcf).

## Summary and conclusion
According to found variants in provided data, the [CARD11](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000198286;r=7:2906142-3044228;tl=B5j5OQNfJqj0ZkRK-10087084) gene in the donor of the DNA is highly variant. Majority of the found variants were located in non-coding parts of the gene. 3 variants were located in the protein-coding parts of the gene, out of which one was stop mutation.

Stop mutation is affecting trancript [ENST00000396946.9](https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000198286;r=7:2906142-3043867;t=ENST00000396946;tl=B5j5OQNfJqj0ZkRK-10087084), in the way, that the way, that transcript is ending prematurely (only 21 out of 25 exons will be transcribed).

This is likely resulting in the disfunction of the transcript.



