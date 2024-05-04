# Cancer Genomics Data Analysis

First, we will, again, prepare virtual environment, where we will work.

### Create environment and get the data

```bash
virtualenv exercise_1
source exercise_1/bin/activate
sudo apt install -y bwa
sudo apt install -y bcftools
sudo apt install -y samtools
sudo apt install -y r-base-core
sudo apt install -y delly
sudo apt install -y tabix

```
Next, we will get our data. Data labeled as *tu* are reads from the tumor sample. Data labeled as *wt* are reads from the germline sample.

```bash
wget https://gear-genomics.embl.de/data/.exercise/tu.r1.fq.gz
wget https://gear-genomics.embl.de/data/.exercise/tu.r2.fq.gz
wget https://gear-genomics.embl.de/data/.exercise/wt.r1.fq.gz
wget https://gear-genomics.embl.de/data/.exercise/wt.r2.fq.gz
```

### Alignment using bwa-mem
Same as in Exercise_1, we will map the reads to the reference genome. As reference genome, we will use
human genome from [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz).

We will create separate alignment to reference genome using tumor and germline samples.
Similarly to previous assignment, we will need to first create index file from reference, before we could start with mapping. 

```bash
# get the reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
# create index files from reference
bwa index hg19.fa.gz

# map the germline sample, using 4 threads for faster computation
bwa mem -t 4 hg19.fa.gz wt.r1.fq.gz wt.r2.fq.gz > alignment_germline.bam

# map the tumor sample, using 4 threads for faster computation
bwa mem -t 4 hg19.fa.gz tu.r1.fq.gz tu.r2.fq.gz > alignment_tumor.bam
```

### Variant calling using delly
Now, we will do variant calling using structural-variant caller delly.
Delly will need .fai index file of the reference genome, so we will generate this first.

Also we will need index files for both alignments - for tumor and germline samples.

```bash
gzip -d hg19.fa.gz
samtools faidx hg19.fa

samtools sort -@ 4 alignment_germline.bam -o alignment_germline_sorted.bam
samtools sort -@ 4 alignment_tumor.bam -o alignment_tumor_sorted.bam
samtools index alignment_tumor_sorted.bam
samtools index alignment_germline_sorted.bam

delly call -q 20 -g hg19.fa -o variants_output.bcf alignment_tumor_sorted.bam alignment_germline_sorted.bam
```
Next, we need to do filtering somatic variants in the output file. For this, we will first generate tsv file `samples.tsv` with following content:
```
alignment_tumor_sorted tumor
alignment_germline_sorted control

```
Then we can continue to sort:
```bash
delly filter -p -f somatic -o somatic.bcf -a 0.25 -s samples.tsv variants_output.bcf
```
