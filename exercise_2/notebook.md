# Cancer Genomics Data Analysis

First, we will, again, prepare virtual environment, where we will work.

### Create environment and get the data

```bash
virtualenv exercise_1
source exercise_1/bin/activate
sudo apt install -y bwa
sudo apt install -y bcftools
sudo apt install -y samtools
sudo apt install r-base-core

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

# map the germline sample
bwa mem hg19.fa.gz wt.r1.gz wt.r2.gz > alignment_germline.bam

# map the tumor sample
bwa mem hg19.fa.gz tu.r1.gz tu.r2.gz > alignment_tumorls.bam
```
