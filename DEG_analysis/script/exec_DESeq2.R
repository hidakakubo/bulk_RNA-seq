# script to perform differential gene expression analysis using DESeq2 package

# clean up R workspace
rm(list=ls(all=TRUE))

# load libraries
library(DESeq2)
library(tidyverse)

# cd
path <- "/Users/hoge/Desktop/bulk_RNA-seq/DEG_analysis"
setwd(path)

# make dir for outputs
dir.create("DESeq2_outputs")

# cd to data dir
setwd(str_interp("${path}/data"))



# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv(str_interp("gene_count_matrix.csv"), check.names=FALSE, row.names = 1)

# read in sample info
colData <- read.csv(str_interp("metadata.csv"), row.names = 1)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))



# Step 2: construct a DESeqDataSet object ----------
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ genotype)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$genotype <- relevel(dds$genotype, ref = "WT")

# NOTE: collapse technical replicates



# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res

# sort by padj
res_ordered <- res[order(res$padj),]

res_ordered

# to csv
write.csv(as.data.frame(res_ordered), 
          file=str_interp("${path}/DESeq2_outputs/result.csv"))



# Explore Results ----------------
# summary
summary(res)

# summary(p-value < 0.05) to txt
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

sink(str_interp("${path}/DESeq2_outputs/summary.txt"))
summary(res0.05)
sink()

# normalized counts to csv
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), 
          file=str_interp("${path}/DESeq2_outputs/normalized_counts.csv"))

# MA plot
png(str_interp("${path}/DESeq2_outputs/MAplot.png"), width = 2000, height =1600, res = 400)
plotMA(res)
dev.off()



# Option Results ----------------
# construct a "log-transformed" DESeqDataSet object
ntd <- normTransform(dds)
ntd

# PCA
png(str_interp("${path}/DESeq2_outputs/PCAplot.png"), width = 2000, height =1600, res = 400)
plotPCA(ntd, intgroup=c("genotype"))
dev.off()