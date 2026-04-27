# Analysis reproduction order

All scripts should be run from the repository root directory.

## Prerequisites
- salmon.merged.gene_counts.tsv in results/salmon/ (output from nf-core/rnaseq pipeline)
- Klaeger matrix at /path/to/klaeger_matrix.csv (from companion selectivity repo)

## Step 1 — Differential expression
Rscript analysis/deseq2_analysis.R

## Step 2 — Pathway activity inference
Rscript analysis/progeny_analysis.R

## Step 3 — Pathway enrichment
Rscript analysis/fgsea_analysis.R

## Step 4 — Binding vs transcription
python3 analysis/binding_vs_transcription.py

## Step 5 — Kinase-substrate analysis
python3 analysis/kinase_substrate_analysis.py

## Step 6 — Expression variability
python3 analysis/expression_variability.py

## Step 7 — Transcriptional specificity
python3 analysis/transcriptional_specificity/specificity_analysis.py
