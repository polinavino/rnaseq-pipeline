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

## Dependencies

The binding affinity analyses (`binding_vs_transcription.py`,
`kinase_substrate_analysis.py`, `transcriptional_specificity/specificity_analysis.py`)
require the Klaeger et al. chemoproteomic dataset, available from the companion
repository at https://github.com/polinavino/kinase-selectivity-definitions
(file: `klaeger_matrix.csv`). Set the path via environment variable:

```bash
export KLAEGER_PATH=/path/to/klaeger_matrix.csv
```

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
and scoring methods, so they can disagree even on the same data. This is a
distinct phenomenon from the definitional instability analyzed in the companion
selectivity paper — here the disagreement arises from different statistical
frameworks applied to the same question, rather than different operationalizations
of the same concept from the same data. Both are instances of the broader problem
that biological conclusions can be sensitive to methodological choices, motivating
the kind of formal analysis the selectivity paper pursues.

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

## Binding affinity vs transcriptional change

**Script:** `analysis/binding_vs_transcription.py`  
**Data:** Klaeger et al. chemoproteomic dataset (343 kinases) + DESeq2 results

For each of the 343 kinases in imatinib's binding profile (Klaeger dataset),
we asked whether binding affinity predicts transcriptional change of that
kinase's own gene, and separately whether it predicts the mean transcriptional
change of that kinase's known substrates (PhosphoSitePlus via Harmonizome).

### Results

**Kinase gene expression:** Pearson r = -0.077, p = 0.19 — no correlation.
Expected, since kinase inhibitors act post-translationally on the protein, not
on gene transcription.

**Substrate expression (kinase-substrate analysis):** Spearman r = -0.050,
p = 0.84 for mean substrate log2FC; r = -0.448, p = 0.054 for mean |log2FC|.
Neither is significant. Only 19 of 50 kinases had sufficient substrate
annotation, and most clustered at pKd = 5.0 (the Klaeger detection threshold),
severely limiting the dynamic range of the test.

### Interpretation

These are null results at the gene and substrate level. They confirm that
binding-based selectivity metrics — whatever definition is used — measure
protein-level target engagement, not transcriptional specificity. The
transcriptional response to kinase inhibition is mediated through signaling
cascades and transcription factors, not through direct regulation of target
gene expression. The appropriate level of analysis is pathway-level, as shown
by the PROGENy and fgsea results above.

## Next steps

* Cross-compound comparison — find RNA-seq data for dasatinib (broader
  selectivity profile than imatinib) and test whether transcriptional breadth
  scales with binding promiscuity
* Formal quantification of transcriptional specificity — develop a
  transcriptional analog of the binding selectivity metrics from the companion
  paper and compare them
* Extension to other CML drugs — ponatinib, bosutinib, nilotinib all have
  Klaeger binding profiles and published RNA-seq data exists for some

## Transcriptional specificity analysis

**Script:** `analysis/transcriptional_specificity/specificity_analysis.py`  
**Full writeup:** `analysis/transcriptional_specificity/README.md`

This analysis applies the binding selectivity desiderata from the companion
paper to PROGENy pathway z-scores, testing whether the same definitional
instabilities appear when selectivity-type measures are computed over
transcriptional rather than binding data.

### Scores

| Metric | Value | Note |
|--------|-------|------|
| S-score analog | 0.571 | 8/14 pathways with \|z\| > 2 |
| Entropy analog | 3.054 bits | Moderate specificity |
| Gini analog | 0.552 | Moderate concentration |
| Ratio analog | 1.288 | MAPK 29% stronger than PI3K |

### Desiderata results

D1 violated by S-score (threshold-sensitive), D3 violated by Gini (decreases
when weak pathway added — opposite to entropy), D3 violated by ratio
(insensitive to off-target pathways). D4 passes for all definitions.

### Significance

The same failure modes identified in binding-based selectivity metrics appear
in transcriptional pathway data. This suggests the desiderata violations are
properties of the mathematical definitions, not artifacts of binding
measurement. The companion paper's framework generalizes beyond binding
selectivity to any domain where specificity is measured over a response profile.

## Expression variability analysis

**Script:** `analysis/expression_variability.py`  
**Motivated by:** Kaern et al. — noise-driven drug resistance framework

This analysis asks whether expression variability across biological replicates
predicts the transcriptional response to imatinib, connecting our bulk RNA-seq
results to the theoretical framework developed by Kaern and colleagues on
stochastic gene expression and non-genetic drug resistance.

### Background

Kaern's lab has shown mathematically and experimentally that genetically
identical cells can develop drug resistance through stochastic variability in
gene expression — cells that happen to express a resistance gene at high levels
when drug is applied survive, without any genetic mutation. This predicts that
genes with high expression variability in untreated cells should be enriched
among genes that respond to drug treatment, since variable genes are more likely
to have subpopulations already primed for survival or death.

### Results

**Q1 — Does control variability predict transcriptional response?**

Spearman r = 0.231 (p ≈ 0) between control coefficient of variation (CV) and
|log2FC| in imatinib treatment. There is a weak but highly significant positive
correlation — more variable genes tend to show larger transcriptional responses.
This is consistent with Kaern's framework: noisy genes are more likely to be
selected for or against by drug treatment. The effect is modest, suggesting
noise is one factor among many.

**Q2 — Does imatinib reduce expression variability?**

Mean CV increases from 0.144 (control) to 0.165 (imatinib), a statistically
significant increase (paired t-test p = 4e-75). This is the opposite of what
a simple clonal selection model would predict — if imatinib were killing the
most variable cells and leaving a uniform surviving population, CV should
decrease. Instead, imatinib appears to increase transcriptional heterogeneity,
consistent with cells being pushed toward erythroid differentiation at different
rates. This aligns with the dominant fgsea finding (HEME_METABOLISM NES = 3.52)
— differentiation is occurring but heterogeneously across the population.

**Q3 — Are high-variability genes enriched among DE genes?**

High-CV genes (top 10% by control CV) are actually *depleted* among significant
DEGs — 47.0% vs 61.8% overall. This is a statistical artifact: high CV genes
have noisy expression that makes consistent DE calling difficult. But it also
has a biological interpretation: the most variable genes may respond to imatinib
inconsistently across replicates, with some cells responding strongly and others
not. Bulk RNA-seq, which averages across millions of cells, cannot resolve this.

### Connection to Kaern's framework

The results here are exploratory and should be interpreted cautiously given the
limitations of bulk RNA-seq for testing noise-driven resistance predictions.

The increase in expression variability under imatinib (Q2) is consistent with
heterogeneous differentiation — cells transitioning toward erythroid fate at
different rates, as seen in the fgsea HEME_METABOLISM result. However, the
effect size is small (+0.022 mean CV change) and bulk RNA-seq cannot distinguish
true cell-to-cell variability from technical replicate noise. The weak positive
correlation between control CV and |log2FC| (r=0.231) is real but modest and
not surprising. The depletion of high-CV genes among DEGs is a statistical
artifact of noisy expression being harder to call as significantly DE.

These analyses motivate single-cell RNA-seq follow-up, which would be needed to
properly test Kaern's predictions about subpopulation-level responses to
imatinib. The SRP562191 dataset (dasatinib-treated K562 Multiome) offers a
partial opportunity, though with only 2 replicates.

### Limitation and next steps

This analysis uses replicate-to-replicate variability as a proxy for
cell-to-cell variability. Each replicate is an average of millions of cells,
so our CV reflects biological variation between experimental replicates plus
technical noise, not true single-cell heterogeneity. Single-cell RNA-seq of
imatinib-treated K562 cells would be needed to directly test Kaern's predictions
about subpopulation-level responses. The SRP562191 dataset (dasatinib-treated
K562 Multiome, 2 replicates) offers a partial opportunity for comparison,
though with fewer replicates than ideal.

## Nilotinib RNA-seq analysis

**Dataset:** SRP594028 — K562 Mock+DMSO (control) vs Mock+Nilotinib, 3 replicates each  
**Pipeline:** Same nf-core/rnaseq configuration as imatinib (GCP, Salmon pseudoalignment)  
**Script:** `analysis/deseq2_nilotinib.R`, `analysis/progeny_nilotinib.R`

Nilotinib is a second-generation BCR-ABL inhibitor, more potent and more
BCR-ABL-selective than imatinib by binding affinity (Klaeger dataset).
This experiment allows direct comparison of two drugs with different binding
selectivity profiles but the same nominal target in the same cell line.

### DESeq2 results

- Genes tested: 15,451
- Significant DEGs (padj < 0.05): 9,478
- Upregulated: 1,881 | Downregulated: 1,833

Known BCR-ABL response genes all show the expected direction. Erythroid
differentiation markers (HBZ log2FC=+5.55, ALAS2=+4.63, HBG2=+3.00) are
more strongly upregulated by nilotinib than imatinib, consistent with
nilotinib's greater BCR-ABL potency.

## Imatinib vs Nilotinib comparison

**Script:** `analysis/compare_imatinib_nilotinib.py`  
**Output:** `analysis/imatinib_vs_nilotinib.png`, `analysis/specificity_comparison.csv`

### PROGENy pathway comparison

| Pathway | Imatinib z | Nilotinib z | Note |
|---------|-----------|------------|------|
| MAPK | -18.84 | -11.44 | Imatinib suppresses more strongly |
| PI3K | -14.62 | -5.02 | Imatinib suppresses more strongly |
| EGFR | -13.44 | -4.19 | Imatinib suppresses more strongly |
| Hypoxia | +4.30 | -4.23 | Opposite directions |
| JAK-STAT | +0.48 | +3.85 | Nilotinib activates more strongly |
| p53 | +2.46 | -0.07 | Opposite directions |
| TGFb | -5.21 | -5.27 | Nearly identical |

### Transcriptional specificity comparison

| Metric | Imatinib | Nilotinib | Agreement |
|--------|----------|-----------|-----------|
| S-score | 0.571 | 0.786 | Nilotinib less specific |
| Entropy | 3.054 | 3.456 | Nilotinib less specific |
| Gini | 0.552 | 0.371 | Nilotinib less specific |
| Ratio | 1.288 | 1.558 | Nilotinib more specific |

Three of four definitions agree: nilotinib produces a broader transcriptional
response than imatinib. Only the ratio definition disagrees — the same D3
violation pattern identified in the companion selectivity paper.

### Key finding

Nilotinib is more potent and more BCR-ABL-selective than imatinib at the
binding level (Klaeger dataset), yet produces a broader transcriptional
response by three of four specificity definitions. MAPK and PI3K suppression
are weaker for nilotinib despite stronger BCR-ABL inhibition. Hypoxia and
p53 pathways show opposite activation directions between the two drugs.

This is direct empirical evidence that binding selectivity does not predict
transcriptional specificity. A more BCR-ABL-selective compound does not
necessarily produce a more transcriptionally specific response. The
relationship between binding profile and transcriptional footprint is
mediated by signaling network architecture in ways not captured by any
binding-based selectivity definition.
