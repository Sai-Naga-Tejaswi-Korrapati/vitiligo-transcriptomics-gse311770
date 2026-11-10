# =====================================================================
# Project: Transcriptomic Profiling of Vitiligo Skin (GSE311770)
# Author: M.Sc. Bioinformatics Portfolio
# Objective: Identify Differentially Expressed Genes (DEGs) between 
#            Lesional and Control skin samples using DESeq2.
# =====================================================================

# 1. Load Required Libraries
library(DESeq2)
library(org.Hs.eg.db)
library(EnhancedVolcano)

# 2. Load Count Matrix
# (Ensure your working directory is set to where the counts file is saved)
#  Pop up a window to select your counts file (choose the .tsv file)
file_path <- file.choose() # Update name if needed
counts <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")

# Convert to a clean numeric matrix
counts_num <- as.matrix(counts)
class(counts_num) <- "numeric"

# 3. Construct Metadata Dynamically from Sample Names
sample_names <- colnames(counts_num)
condition <- ifelse(grepl("lesioneel", sample_names) & !grepl("perilesioneel", sample_names), "Lesional",
             ifelse(grepl("perilesioneel", sample_names), "Perilesional",
             ifelse(grepl("gezond", sample_names), "NonLesional", "Control")))

meta <- data.frame(
  row.names = sample_names,
  condition = as.factor(condition)
)

# Verify group distributions
table(meta$condition)

# 4. Build DESeqDataSet and Filter Low-Count Genes
dds <- DESeqDataSetFromMatrix(countData = round(counts_num), 
                              colData = meta, 
                              design = ~ condition)

# Keep genes with at least 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 5. Run Differential Expression Analysis
dds <- DESeq(dds)

# Extract results for Lesional vs. Control
res <- results(dds, contrast = c("condition", "Lesional", "Control"))

# 6. Annotate ENSEMBL IDs with Gene Symbols
res_df <- as.data.frame(res)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# 7. Filter Significant DEGs (Padj < 0.05 & |log2FC| > 1)
sig_res <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
sig_res <- sig_res[order(sig_res$padj), ]

# Export significant results to CSV
write.csv(sig_res, file = "Vitiligo_Lesional_vs_Control_SigDEGs.csv", row.names = TRUE)

# 8. Generate Publication-Quality Volcano Plot
pdf("Vitiligo_Volcano_Plot.pdf", width = 8, height = 7)
EnhancedVolcano(res_df,
                lab = res_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 4.0,
                title = 'Vitiligo Lesional vs. Control',
                subtitle = 'GSE311770 RNA-Seq Analysis',
                legendPosition = 'right')
dev.off()
cat("Analysis complete! Significant DEGs saved and volcano plot generated successfully.\n")

# 1. Print the top 10 significant genes directly to your console
head(sig_res, 10)

# 2. Open a spreadsheet-style viewer window for all your significant genes
View(sig_res)

# 3. Pop up the Volcano plot immediately in your RStudio Plots window
library(EnhancedVolcano)
EnhancedVolcano(res_df,
                lab = res_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 4.0,
                title = 'Vitiligo Lesional vs. Control',
                subtitle = 'GSE311770 RNA-Seq Analysis')