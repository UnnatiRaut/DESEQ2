library(DESeq2)
library(ggplot2)
counts <- read.csv("IPF_stages//GSE213001//Early_advanced_124_M.csv", row.names=1)
coldata <- read.csv("IPF_stages//GSE213001//metadata_M.csv", row.names=1)

colnames(coldata)

coldata$Patient <- as.factor(coldata$Patient)
coldata$Condition <- as.factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Patient + Condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Advanced", "Early"))
res <- res[order(res$padj), ]
sig_genes <- subset(res, padj < 0.05)
write.csv(as.data.frame(res), "DESeq2_results_paired.csv")


vsd <- vst(dds, blind=FALSE)  # Use variance-stabilizing transformation
plotPCA(vsd, intgroup = c("Condition"))
pca_data <- plotPCA(vsd, intgroup = c("Condition", "Patient"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar")) # Get percentage variance
shape_values <- c(16, 17, 18, 19, 15, 7, 8, 9, 10, 11, 12, 13, 14, 21, 22, 23, 24, 25, 4, 3)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Patient)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA of Lung Samples") +
  theme_minimal() +
  scale_color_manual(values = c("Apex" = "blue", "Base" = "red")) +
  scale_shape_manual(values = shape_values[1:length(unique(pca_data$Patient))])  # Assign shapes dynamically

# Get normalized counts from DESeq2
normalized_counts <- counts(dds, normalized = TRUE)

# Convert to a data frame
normalized_counts_df <- as.data.frame(normalized_counts)

# Save the normalized counts to a CSV file
write.csv(normalized_counts_df, "DESeq2_normalized_counts.csv")
