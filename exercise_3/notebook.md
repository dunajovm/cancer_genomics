# Working with count matrices 

### Get the data
Firs, we will get the data

```bash
wget https://gear-genomics.embl.de/data/.slides/sample.counts
wget https://gear-genomics.embl.de/data/.slides/sample.info
wget https://gear-genomics.embl.de/data/.slides/template.R

```

### Explore the template
Lets open the template in the Rstudio program.

```R
library(DESeq2)
library(pheatmap)

# Load count matrix
x = read.table("sample.counts", row.names=1, header=T, sep=",")
s = read.table("sample.info", header=T, row.names=1, colClasses=c("character", "factor"))

# Create DESeq2 object
dds = DESeqDataSetFromMatrix(countData = x, colData = s, design = ~ condition)
```

### Differential expression analysis
Run a differential expression analysis (Tumour vs. Normal) using a log-fold change threshold of 1.

```R

```

### Generating PCA

### Generating MA-plots

### Exporting into CSV
