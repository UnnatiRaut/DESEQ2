library(DESeq2)
library(tidyverse)
library(airway)
library(writexl)

sampleData <- read.csv("sample FC\\sampleData.csv",row.names = 1)
Degs_data <- read.csv("sample FC\\DEG data_1.csv",row.names = 1)
head(sampleData)
head(Degs_data)
print(all(colnames(Degs_data) %in% rownames(sampleData)))
all(colnames(Degs_data) == rownames(sampleData))
Degs_data <- na.omit(Degs_data)

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              sampleData = sampleData,
                              design = ~ Dexamethasome)
dds
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

dds

dds$Dexamethasome <- relevel(dds$Dexamethasome, ref = "control")

dds <- DESeq(dds)
res <- results(dds)

res
summary(res)

res0.01 <- results(dds, alpha = 0.1)
summary(res0.01)

# contrasts
resultsNames(dds)

plotMA(res)


res <- results(dds)
res_df <- as.data.frame(res)
write_xlsx(res_df, "DESeq2_results.xlsx")
