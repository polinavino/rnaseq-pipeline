library(progeny)
library(dplyr)
library(tibble)
library(readr)

res_df <- read_csv("analysis/deseq2_nilotinib_results.csv")

gene_stats <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(stat = sign(log2FoldChange) * -log10(padj + 1e-300)) %>%
    dplyr::select(gene_id, stat) %>%
    deframe()

cat("Genes for PROGENy:", length(gene_stats), "\n")

pathway_scores <- progeny(
    expr     = matrix(gene_stats, ncol=1,
                      dimnames=list(names(gene_stats), "nilotinib_vs_control")),
    scale    = FALSE,
    organism = "Human",
    top      = 500,
    perm     = 1000,
    z_scores = TRUE
)

scores_df <- as.data.frame(t(pathway_scores)) %>%
    rownames_to_column("pathway") %>%
    rename(z_score_nilotinib = nilotinib_vs_control) %>%
    arrange(z_score_nilotinib)

cat("\nNilotinib PROGENy z-scores:\n")
print(scores_df)

write_csv(scores_df, "analysis/progeny_nilotinib_scores.csv")
