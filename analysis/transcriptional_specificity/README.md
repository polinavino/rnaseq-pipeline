# Transcriptional Specificity Analysis

This analysis applies the binding selectivity desiderata from the companion
paper (ChemRxiv: https://doi.org/10.26434/chemrxiv.15001618/v1) to
transcriptional pathway data, asking whether the same definitional instabilities
appear when selectivity-like measures are applied to a different data type.

## Motivation

The companion paper identified four desiderata (D1-D4) that any well-formed
kinase inhibitor selectivity measure should satisfy, and showed that no existing
binding-based definition satisfies all four simultaneously. Here we ask: are
these failures specific to binding data, or do they reflect a more general
problem with selectivity-type measurement?

## Method

We treat imatinib's PROGENy pathway activity profile (14 pathways, z-scores
from 1000 permutations) as an analog of a binding profile, and compute four
transcriptional specificity scores:

- **S-score analog**: fraction of pathways with |z| > threshold
- **Entropy analog**: Shannon entropy of normalized |z| distribution
- **Gini analog**: Gini coefficient of |z| distribution
- **Ratio analog**: top pathway |z| / second pathway |z|

We then test each score against D1-D4.

## Results

### Transcriptional specificity scores

| Metric | Value | Interpretation |
|--------|-------|----------------|
| S-score | 0.571 | 8/14 pathways strongly perturbed (|z| > 2) |
| Entropy | 3.054 bits | Moderate (max = 3.81 bits for 14 pathways) |
| Gini | 0.552 | Moderate concentration |
| Ratio | 1.288 | MAPK only 29% stronger than PI3K |

### Desiderata results

**D1 (reliability threshold):** The S-score is highly threshold-sensitive
(Panel E). At |z|=0.5, 93% of pathways appear perturbed. At |z|=4.5, only
43% do. The choice of threshold has no principled basis — the same violation
as in the binding setting.

**D2 (bounded gap sensitivity):** The ratio captures a 29% premium for MAPK
over PI3K. Entropy sees this as 0.37 bits. Both are modest, and both tools
agree the gap is small — D2 is approximately satisfied in this case.

**D3 (distributional consistency):** Adding a weak off-target pathway (|z|=0.5):
- Entropy increases (correctly registers less specificity) ✓
- Gini *decreases* (incorrectly registers more specificity) ✗ — D3 violation
- Ratio is unchanged (insensitive to off-target pathways below top two) ✗ — D3 violation
- S-score decreases due to denominator growth — D3 violation

Entropy and Gini give opposite responses to the same perturbation (Panel C),
reproducing the binding paper's finding that these two definitions measure
fundamentally different properties.

**D4 (monotonicity):** All four definitions pass. Strengthening the weakest
pathway reduces specificity under all metrics.

### Key finding

The same definitional instabilities identified in binding-based selectivity
metrics appear when the framework is applied to transcriptional pathway data.
Gini and entropy disagree on D3; ratio is insensitive to distributional changes;
S-score is threshold-dependent. These failures are not specific to binding data
— they reflect properties of the mathematical definitions themselves, independent
of the underlying data type.

This supports the generalization of the companion paper's desiderata framework
beyond binding-based selectivity to any domain where a "specificity" measure
is computed over a profile of responses across a set of targets or pathways.
