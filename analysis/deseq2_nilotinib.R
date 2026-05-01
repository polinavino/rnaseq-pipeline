library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(readr)

counts <- read_tsv("results/nilotinib/salmon/salmon.merged.gene_counts.tsv") %>%
    column_to_rownames("gene_id") %>%
    select(-gene_name) %>%
    round() %>%
    as.matrix()

coldata <- data.frame(
    sample    = colnames(counts),
    condition = factor(c("control", "control", "control",
                         "nilotinib", "nilotinib", "nilotinib"),
                       levels = c("control", "nilotinib")),
    row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","nilotinib","control"),
               alpha=0.05, pAdjustMethod="BH")

res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

cat("Significant genes (padj<0.05):", sum(res_df$padj < 0.05, na.rm=TRUE), "\n")
cat("Upregulated:", sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm=TRUE), "\n")
cat("Downregulated:", sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm=TRUE), "\n")

write_csv(res_df, "analysis/deseq2_nilotinib_results.csv")

vsd <- vst(dds, blind=FALSE)

p_pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE) %>%
    {ggplot(., aes(PC1, PC2, color=condition, label=name)) +
     geom_point(size=4) + geom_text(vjust=-0.8, size=3) +
     scale_color_manual(values=c("control"="steelblue","nilotinib"="tomato")) +
     theme_bw() + ggtitle("PCA: Control vs Nilotinib-treated K562 cells")}
ggsave("analysis/pca_nilotinib.png", p_pca, width=7, height=5, dpi=150)

known_genes <- c("HBZ","ALAS2","HBG1","HBG2","MYC","BCL2","CCND1","CDK4","E2F1","MCM2")
known_results <- res_df %>% filter(gene_id %in% known_genes) %>%
    select(gene_id, log2FoldChange, padj)
cat("\nKnown BCR-ABL response genes:\n")
print(known_results)

write_csv(known_results, "analysis/nilotinib_known_gene_results.csv")
cat("\nDone.\n")
