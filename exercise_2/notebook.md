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

```bash
# get the reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```
