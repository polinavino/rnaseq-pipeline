# Bulk RNA-seq Pipeline: Kinase Inhibitor Treatment Response

Transcriptional response of K562 CML cells to imatinib treatment, analyzed
using nf-core/rnaseq on Google Cloud Platform and DESeq2.

## Biological question

How does imatinib — a BCR-ABL kinase inhibitor — reshape the transcriptome
of K562 chronic myeloid leukemia cells? Which pathways are activated or
suppressed, and are the transcriptional changes consistent with imatinib's
known kinase binding selectivity profile?

This project connects to a companion study on kinase inhibitor selectivity
definitions (ChemRxiv: https://doi.org/10.26434/chemrxiv.15001618/v1), which
characterizes imatinib's binding profile across 343 kinases using the Klaeger
chemoproteomic dataset.

## Dataset

**GEO accession:** SRP528723  
**Cell line:** K562 (human CML, BCR-ABL+)  
**Treatment:** Imatinib (1 µM) vs. DMSO control  
**Replicates:** 3 per condition  
**Sequencing:** Paired-end, Illumina NovaSeq 6000, ~21M read pairs/sample

| Sample | SRR accession |
|--------|--------------|
| control_rep1 | SRR30403510 |
| control_rep2 | SRR30403509 |
| control_rep3 | SRR30403508 |
| imatinib_rep1 | SRR30403504 |
| imatinib_rep2 | SRR30403503 |
| imatinib_rep3 | SRR30403502 |

## Pipeline

**Workflow:** nf-core/rnaseq v3.24.0  
**Executor:** Google Cloud Batch (us-central1)  
**Reference genome:** GRCh38 (iGenomes)  
**Quantification:** Salmon (pseudo-alignment, --skip_alignment)

### Key pipeline steps
- FastQC — raw read quality control
- TrimGalore — adapter trimming
- Salmon — transcript quantification
- tximeta/tximport — count aggregation
- MultiQC — QC report aggregation

### Infrastructure
- Google Cloud Storage bucket for FASTQs, work directory, and results
- Google Batch for job execution
- Nextflow v25.10.4 for workflow orchestration

## Differential expression analysis

**Script:** `analysis/deseq2_analysis.R`  
**Tool:** DESeq2  
**Comparison:** imatinib vs. control  
**Significance threshold:** padj < 0.05, |log2FC| > 1

### Results summary
- Genes tested: 13,480 (after low-count filtering)
- Significant DEGs: 8,319
- Upregulated: 1,527
- Downregulated: 991

### Known imatinib response genes

| Gene | log2FC | padj | Interpretation |
|------|--------|------|----------------|
| HBZ | +2.99 | 9e-96 | Erythroid differentiation ↑ |
| ALAS2 | +3.46 | 5e-35 | Erythroid differentiation ↑ |
| HBG2 | +1.52 | 5e-35 | Erythroid differentiation ↑ |
| HBG1 | +1.39 | 2e-29 | Erythroid differentiation ↑ |
| BCL2 | -1.62 | 2e-09 | Pro-survival suppressed |
| MYC | -1.51 | 3e-08 | Proliferation suppressed |
| CCND1 | -2.33 | 4e-05 | Cell cycle arrest |

Results are consistent with BCR-ABL inhibition driving K562 cells toward
erythroid differentiation and away from proliferation.

## Outputs

| File | Description |
|------|-------------|
| `analysis/deseq2_results.csv` | Full DESeq2 results table |
| `analysis/known_gene_results.csv` | Known imatinib response genes |
| `analysis/pca.png` | PCA plot |
| `analysis/volcano.png` | Volcano plot |
| `analysis/heatmap_top50.png` | Heatmap of top 50 DEGs |

## Reproduction

### Requirements
- Nextflow >= 25.0
- Google Cloud SDK
- R >= 4.0 with DESeq2, ggplot2, pheatmap, readr, dplyr, tibble

### GCP setup
```bash
gcloud auth login
gcloud config set project YOUR_PROJECT_ID
gsutil mb -l us-central1 gs://YOUR_BUCKET
```

### Run pipeline
```bash
export GOOGLE_APPLICATION_CREDENTIALS=~/nextflow-sa-key.json

nextflow run nf-core/rnaseq \
  -profile googlebatch \
  -c custom.config \
  --input gs://YOUR_BUCKET/samplesheet.csv \
  --outdir gs://YOUR_BUCKET/results \
  --genome GRCh38 \
  --skip_linting \
  --skip_alignment \
  --pseudo_aligner salmon \
  --project_id YOUR_PROJECT_ID \
  --workdir_bucket gs://YOUR_BUCKET/work \
  --workers_service_account YOUR_SA@YOUR_PROJECT.iam.gserviceaccount.com \
  --use_spot false \
  --boot_disk "50 GB" \
  -queue-size 2
```

### Run DESeq2 analysis
```bash
Rscript analysis/deseq2_analysis.R
```

## Next steps

- Kinase activity inference using PROGENy
- Pathway enrichment analysis (fgsea)
- Cross-reference inferred kinase activities with imatinib binding profile
  from the Klaeger dataset
- Comparison with dasatinib transcriptional response (broader selectivity
  profile, predicted broader transcriptional footprint)


## Kinase activity inference (PROGENy)

**Script:** `analysis/progeny_analysis.R`  
**Tool:** PROGENy (Pathway RespOnsive GENes)  
**Permutations:** 1000  
**Top genes per pathway:** 500

PROGENy infers upstream pathway activity from the DESeq2 gene-level statistics,
using a database of pathway-responsive genes derived from perturbation experiments.

### Pathway z-scores

| Pathway | Z-score | Interpretation |
|---------|---------|----------------|
| MAPK | -18.84 | Strongly suppressed — primary BCR-ABL downstream pathway |
| PI3K | -14.62 | Strongly suppressed — primary BCR-ABL downstream pathway |
| EGFR | -13.44 | Strongly suppressed — imatinib does not bind EGFR |
| Estrogen | -7.97 | Suppressed |
| TGFb | -5.21 | Suppressed |
| VEGF | -1.55 | Weakly suppressed |
| WNT | -1.28 | Weakly suppressed |
| NFkB | -1.13 | Weakly suppressed |
| Androgen | -1.01 | Weakly suppressed |
| TNFa | +0.08 | Unchanged |
| JAK-STAT | +0.48 | Unchanged (expected suppression not observed) |
| p53 | +2.46 | Activated — BCR-ABL inhibition relieves p53 suppression |
| Trail | +3.69 | Activated — apoptotic signaling |
| Hypoxia | +4.30 | Strongly activated |

### Connection to binding selectivity

Imatinib's primary targets in the Klaeger chemoproteomic dataset are BCR-ABL,
c-KIT, and PDGFR. MAPK and PI3K suppression (z = -18.8 and -14.6) are
consistent with BCR-ABL inhibition. However, the EGFR pathway is equally
suppressed (z = -13.4) despite imatinib having no meaningful binding affinity
for EGFR in the Klaeger dataset.

This illustrates a key finding from the companion selectivity paper: binding
selectivity profiles do not map cleanly onto transcriptional pathway effects.
Pathway crosstalk, shared downstream targets, and compensatory signaling mean
that a compound's transcriptional footprint is broader than its binding profile
would predict — regardless of which selectivity definition is used to
characterize that profile.

JAK-STAT is also notable: despite BCR-ABL being a known activator of JAK-STAT
signaling, the pathway shows near-zero activity change (z = +0.48), suggesting
compensatory activation from other sources in K562 cells.

## Pathway enrichment analysis (fgsea)

**Script:** `analysis/fgsea_analysis.R`  
**Tool:** fgsea with MSigDB Hallmark gene sets  
**Ranked by:** DESeq2 Wald statistic

### Key findings

**Most strongly activated pathways:**

| Pathway | NES | Interpretation |
|---------|-----|----------------|
| HEME_METABOLISM | +3.52 | Erythroid differentiation — dominant imatinib effect in K562 |
| BILE_ACID_METABOLISM | +1.96 | Metabolic reprogramming |
| KRAS_SIGNALING_DN | +1.95 | RAS pathway suppression (downregulated KRAS targets upregulated) |
| HYPOXIA | +1.42 | Metabolic stress response |

**Most strongly suppressed pathways:**

| Pathway | NES | Interpretation |
|---------|-----|----------------|
| MYC_TARGETS_V1 | -2.73 | Proliferation — strongest suppressed signal, BCR-ABL drives MYC |
| MYC_TARGETS_V2 | -2.49 | Proliferation |
| UNFOLDED_PROTEIN_RESPONSE | -2.29 | ER stress pathway suppressed |
| MTORC1_SIGNALING | -2.25 | mTOR pathway — downstream of PI3K/BCR-ABL |
| IL6_JAK_STAT3_SIGNALING | -1.78 | JAK-STAT suppressed (contrast with PROGENy result) |

### Interpretation

The dominant transcriptional effect of imatinib in K562 cells is suppression
of MYC-driven proliferation and activation of erythroid differentiation
(HEME_METABOLISM NES = 3.52). This is consistent with BCR-ABL inhibition
being the primary mechanism — BCR-ABL drives MYC expression, and inhibiting
it collapses the MYC transcriptional program.

The fgsea result for JAK-STAT (NES = -1.78, suppressed) differs from the
PROGENy result (z = +0.48, near zero). The two tools use different gene sets
and scoring methods, illustrating that pathway inference is tool-dependent —
a point relevant to the companion selectivity paper's argument about
definitional instability in kinase biology.

### Connection to binding selectivity

The transcriptional response is highly concentrated: HEME_METABOLISM (NES 3.52)
dominates activated pathways, and MYC_TARGETS (NES -2.73) dominates suppressed
pathways. This concentration is consistent with imatinib's high binding
selectivity for BCR-ABL in the Klaeger dataset — a selective compound produces
a focused transcriptional response.

However, EGFR pathway suppression (PROGENy z = -13.4) in the absence of
EGFR binding demonstrates that binding-based selectivity definitions do not
fully predict pathway-level transcriptional specificity. Pathway crosstalk
and shared downstream targets mean a compound's transcriptional footprint
is broader than its binding profile alone would suggest — a limitation of
any binding-based selectivity framework, regardless of which definition is used.
