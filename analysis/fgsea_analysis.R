library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

# ── 1. Load DESeq2 results ────────────────────────────────────────────────────
res_df <- read_csv("analysis/deseq2_results.csv")

# Ranked gene list by stat (Wald statistic from DESeq2)
ranks <- res_df %>%
    filter(!is.na(stat)) %>%
    dplyr::select(gene_id, stat) %>%
    deframe()

cat("Genes in ranked list:", length(ranks), "\n")

# ── 2. Get gene sets from MSigDB ──────────────────────────────────────────────
# Hallmark gene sets — 50 well-defined biological states
hallmarks <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(f = .$gs_name) %>%
    lapply(function(x) x$gene_symbol)

# KEGG pathways
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(f = .$gs_name) %>%
    lapply(function(x) x$gene_symbol)

# ── 3. Run fgsea ──────────────────────────────────────────────────────────────
set.seed(42)

cat("Running fgsea on Hallmark gene sets...\n")
fgsea_hallmarks <- fgsea(
    pathways = hallmarks,
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500
) %>% arrange(pval)

cat("Running fgsea on KEGG pathways...\n")
fgsea_kegg <- fgsea(
    pathways = kegg,
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500
) %>% arrange(pval)

# ── 4. Save results ───────────────────────────────────────────────────────────
fgsea_hallmarks_out <- fgsea_hallmarks %>%
    dplyr::select(pathway, pval, padj, NES, size) %>%
    arrange(NES)

fgsea_kegg_out <- fgsea_kegg %>%
    dplyr::select(pathway, pval, padj, NES, size) %>%
    arrange(NES)

write_csv(fgsea_hallmarks_out, "analysis/fgsea_hallmarks.csv")
write_csv(fgsea_kegg_out, "analysis/fgsea_kegg.csv")

# ── 5. Print top results ──────────────────────────────────────────────────────
cat("\nTop suppressed Hallmark pathways (imatinib vs control):\n")
print(fgsea_hallmarks_out %>% filter(padj < 0.05) %>% head(10))

cat("\nTop activated Hallmark pathways:\n")
print(fgsea_hallmarks_out %>% filter(padj < 0.05) %>% tail(10))

# ── 6. Plot top Hallmark pathways ─────────────────────────────────────────────
top_pathways <- fgsea_hallmarks_out %>%
    filter(padj < 0.05) %>%
    slice(c(1:10, (n()-9):n())) %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway),
           pathway = factor(pathway, levels = pathway))

p <- ggplot(top_pathways, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                      guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_bw() +
    labs(title = "fgsea Hallmark pathways — Imatinib vs Control (K562)",
         subtitle = "Top 10 suppressed and activated (padj < 0.05)",
         x = "",
         y = "Normalized Enrichment Score (NES)")

ggsave("analysis/fgsea_hallmarks.png", p, width = 10, height = 7, dpi = 150)

cat("\nOutput files:\n")
cat("  analysis/fgsea_hallmarks.csv\n")
cat("  analysis/fgsea_kegg.csv\n")
cat("  analysis/fgsea_hallmarks.png\n")

# ── Fix plot with unique pathway names ───────────────────────────────────────
fgsea_hallmarks_out2 <- read_csv("analysis/fgsea_hallmarks.csv")

top_pathways <- bind_rows(
    fgsea_hallmarks_out2 %>% filter(padj < 0.05) %>% arrange(NES) %>% head(10),
    fgsea_hallmarks_out2 %>% filter(padj < 0.05) %>% arrange(NES) %>% tail(10)
) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway))

p2 <- ggplot(top_pathways, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                      guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_bw() +
    labs(title = "fgsea Hallmark pathways — Imatinib vs Control (K562)",
         subtitle = "Top suppressed and activated (padj < 0.05)",
         x = "",
         y = "Normalized Enrichment Score (NES)")

ggsave("analysis/fgsea_hallmarks.png", p2, width = 10, height = 7, dpi = 150)
cat("Plot saved\n")
