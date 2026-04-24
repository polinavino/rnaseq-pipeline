library(DESeq2)
library(tximeta)
library(tximport)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(readr)

# ── 1. Load data ──────────────────────────────────────────────────────────────
# nf-core/rnaseq outputs a salmon quantification directory per sample
# and a merged gene counts matrix at results/salmon/salmon.merged.gene_counts.tsv

counts_file <- "results/salmon/salmon.merged.gene_counts.tsv"
tpm_file    <- "results/salmon/salmon.merged.gene_tpm.tsv"

counts <- read_tsv(counts_file) %>%
    column_to_rownames("gene_id") %>%
    select(-gene_name) %>%
    round() %>%
    as.matrix()

# ── 2. Sample metadata ────────────────────────────────────────────────────────
coldata <- data.frame(
    sample    = colnames(counts),
    condition = factor(c("control", "control", "control",
                         "imatinib", "imatinib", "imatinib"),
                       levels = c("control", "imatinib")),
    row.names = colnames(counts)
)

stopifnot(all(rownames(coldata) == colnames(counts)))

# ── 3. DESeq2 ─────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,
    design    = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds,
               contrast  = c("condition", "imatinib", "control"),
               alpha     = 0.05,
               pAdjustMethod = "BH")

res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm=TRUE), "\n")
cat("Upregulated:  ", sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm=TRUE), "\n")
cat("Downregulated:", sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm=TRUE), "\n")

write_csv(res_df, "analysis/deseq2_results.csv")

# ── 4. PCA ────────────────────────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.8, size = 3) +
    xlab(paste0("PC1: ", pca_var[1], "% variance")) +
    ylab(paste0("PC2: ", pca_var[2], "% variance")) +
    scale_color_manual(values = c("control" = "steelblue", "imatinib" = "tomato")) +
    theme_bw() +
    ggtitle("PCA: Control vs Imatinib-treated K562 cells")

ggsave("analysis/pca.png", p_pca, width = 7, height = 5, dpi = 150)

# ── 5. Volcano plot ───────────────────────────────────────────────────────────
res_plot <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
        sig = case_when(
            padj < 0.05 & log2FoldChange >  1 ~ "Up",
            padj < 0.05 & log2FoldChange < -1 ~ "Down",
            TRUE ~ "NS"
        )
    )

p_volcano <- ggplot(res_plot, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Up" = "tomato", "Down" = "steelblue", "NS" = "grey70")) +
    theme_bw() +
    labs(title = "Imatinib vs Control — K562 cells",
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value",
         color = "")

ggsave("analysis/volcano.png", p_volcano, width = 7, height = 5, dpi = 150)

# ── 6. Heatmap of top 50 DE genes ─────────────────────────────────────────────
top50 <- res_df %>%
    filter(padj < 0.05) %>%
    slice_max(abs(log2FoldChange), n = 50) %>%
    pull(gene_id)

mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)

annotation_col <- data.frame(
    condition = coldata$condition,
    row.names = coldata$sample
)

ann_colors <- list(condition = c(control = "steelblue", imatinib = "tomato"))

pheatmap(mat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 7,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         filename = "analysis/heatmap_top50.png",
         width = 8, height = 10)

# ── 7. Known imatinib response genes ─────────────────────────────────────────
# These are genes known to respond to BCR-ABL inhibition in K562 cells
known_genes <- c("HBZ", "ALAS2", "HBG1", "HBG2", "MYC", "BCL2",
                 "CCND1", "CDK4", "E2F1", "MCM2")

known_results <- res_df %>%
    filter(gene_id %in% known_genes) %>%
    select(gene_id, log2FoldChange, padj)

cat("\nKnown imatinib response genes:\n")
print(known_results)

write_csv(known_results, "analysis/known_gene_results.csv")

cat("\nAnalysis complete. Output files:\n")
cat("  analysis/deseq2_results.csv\n")
cat("  analysis/pca.png\n")
cat("  analysis/volcano.png\n")
cat("  analysis/heatmap_top50.png\n")
cat("  analysis/known_gene_results.csv\n")
