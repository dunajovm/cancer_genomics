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
dds <- DESeq(dds)
res <- results(dds, lfcThreshold = 1)
res
```

### Generating PCA

```R
vstcounts <- vst(dds, blind=FALSE)
plotPCA(vstcounts, intgroup=c("Status", "CellType"))
```

### Generating MA-plots
The plot visualizes the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales, then plotting these values.

Points will be colored blue if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
```R
plotMA(res, ylim=c(-2,2))
```
```R
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```

### Exporting into CSV

