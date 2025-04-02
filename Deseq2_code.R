library(DESeq2)
library(tidyverse)
library(openxlsx)

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('IPF_ILD//combined//FC_DEGs.csv',row.names = 1)
head(counts_data)


# read in sample info
colData <- read.csv('IPF_ILD//combined//coldata.csv',row.names = 1)
head(colData)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

#dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "control")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

resultsNames(dds)

# Example contrast: treated vs control
res <- results(dds, contrast = c("dexamethasone", "treated", "control"))

# Add gene names (rownames from counts_data) to the results
res$Gene <- rownames(res)

# Save results to Excel
write.xlsx(as.data.frame(res), "GSE199949_DEGs.xlsx")

# MA plot
plotMA(res)