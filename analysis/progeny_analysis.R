library(progeny)
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

# ── 1. Load DESeq2 results ────────────────────────────────────────────────────
res_df <- read_csv("analysis/deseq2_results.csv")

gene_stats <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(stat = sign(log2FoldChange) * -log10(padj + 1e-300)) %>%
    dplyr::select(gene_id, stat) %>%
    deframe()

cat("Genes for PROGENy:", length(gene_stats), "\n")

# ── 2. Run PROGENy with z-scores ──────────────────────────────────────────────
pathway_scores <- progeny(
    expr     = matrix(gene_stats, ncol = 1,
                      dimnames = list(names(gene_stats), "imatinib_vs_control")),
    scale    = FALSE,
    organism = "Human",
    top      = 500,
    perm     = 1000,
    z_scores = TRUE
)

scores_df <- as.data.frame(t(pathway_scores)) %>%
    rownames_to_column("pathway") %>%
    rename(z_score = imatinib_vs_control) %>%
    arrange(z_score)

cat("\nPROGENy pathway z-scores (imatinib vs control):\n")
print(scores_df)

write_csv(scores_df, "analysis/progeny_scores.csv")

# ── 3. Plot ───────────────────────────────────────────────────────────────────
scores_df$pathway <- factor(scores_df$pathway, levels = scores_df$pathway)

p <- ggplot(scores_df, aes(x = pathway, y = z_score, fill = z_score > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                      guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_bw() +
    labs(title = "PROGENy pathway activity — Imatinib vs Control (K562)",
         subtitle = "Z-scores from 1000 permutations",
         x = "",
         y = "Pathway activity z-score")

ggsave("analysis/progeny_pathways.png", p, width = 8, height = 6, dpi = 150)

# ── 4. Interpretation ─────────────────────────────────────────────────────────
cat("\n── Connection to imatinib binding selectivity (Klaeger dataset) ──\n")
cat("Primary target:  BCR-ABL -> MAPK z =", 
    scores_df %>% filter(pathway=="MAPK") %>% pull(z_score) %>% round(2), "\n")
cat("Primary target:  BCR-ABL -> PI3K z =",
    scores_df %>% filter(pathway=="PI3K") %>% pull(z_score) %>% round(2), "\n")
cat("Non-target:      EGFR    -> EGFR z =",
    scores_df %>% filter(pathway=="EGFR") %>% pull(z_score) %>% round(2), 
    "(unexpectedly suppressed)\n")
cat("Non-target:      -       -> JAK-STAT z =",
    scores_df %>% filter(pathway=="JAK.STAT") %>% pull(z_score) %>% round(2),
    "(expected suppression not observed)\n")
