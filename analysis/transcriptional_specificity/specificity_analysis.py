import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr, spearmanr

# ── 1. PROGENy z-scores ───────────────────────────────────────────────────────
# Read from progeny_analysis.R output
progeny_df = pd.read_csv('analysis/progeny_scores.csv')
progeny = pd.Series(
    progeny_df['z_score'].values,
    index=progeny_df['pathway'].values
)

# Use absolute values — we care about magnitude of perturbation, not direction
# (both strongly suppressed and strongly activated = specific response)
abs_z = progeny.abs().sort_values(ascending=False)

print("PROGENy absolute z-scores (ranked):")
print(abs_z.round(3))
print()

# ── 2. Define selectivity metrics ─────────────────────────────────────────────

def s_score(values, threshold=2.0):
    """Fraction of pathways with |z| above threshold.
    Analog of S-score: fraction of kinases above pKd threshold.
    Lower = more specific (fewer pathways perturbed)."""
    return np.sum(values > threshold) / len(values)

def entropy_score(values, beta=0.0):
    """Shannon entropy of normalized |z| distribution.
    Analog of selectivity entropy.
    Lower = more specific (concentrated response)."""
    shifted = np.maximum(values - beta, 0)
    total = np.sum(shifted)
    if total == 0:
        return 0.0
    p = shifted / total
    p = p[p > 0]
    return -np.sum(p * np.log2(p))

def gini_score(values):
    """Gini coefficient of |z| distribution.
    Analog of Gini coefficient.
    Higher = more specific (one pathway dominates)."""
    n = len(values)
    sorted_vals = np.sort(values)
    indices = np.arange(1, n + 1)
    return (2 * np.sum(indices * sorted_vals)) / (n * np.sum(sorted_vals)) - (n + 1) / n

def ratio_score(values):
    """Ratio of top pathway to second pathway |z|.
    Analog of ratio-based selectivity.
    Higher = more specific."""
    sorted_vals = np.sort(values)[::-1]
    if sorted_vals[1] == 0:
        return np.inf
    return sorted_vals[0] / sorted_vals[1]

# ── 3. Compute scores ─────────────────────────────────────────────────────────
vals = abs_z.values

s    = s_score(vals, threshold=2.0)
h    = entropy_score(vals)
g    = gini_score(vals)
r    = ratio_score(vals)

print("── Transcriptional specificity scores (imatinib, K562) ──")
print(f"S-score analog:   {s:.4f}  (fraction of pathways with |z| > 2)")
print(f"Entropy analog:   {h:.4f}  (bits; lower = more specific)")
print(f"Gini analog:      {g:.4f}  (higher = more specific)")
print(f"Ratio analog:     {r:.4f}  (top / second pathway)")
print()

# ── 4. Test desiderata ────────────────────────────────────────────────────────
print("── Desiderata tests ──")

# D1: Reliability threshold
# How many pathways are above the detection threshold (|z| > 1)?
n_active = np.sum(vals > 1.0)
print(f"D1 (reliability): {n_active}/14 pathways have |z| > 1 — score is {'reliable' if n_active >= 3 else 'unreliable'}")

# D2: Bounded gap sensitivity
# Does the ratio definition give unbounded credit to the top pathway?
top1, top2 = sorted(vals)[-1], sorted(vals)[-2]
print(f"D2 (gap sensitivity): top={top1:.2f}, second={top2:.2f}, ratio={top1/top2:.3f}")
print(f"   Ratio captures {((top1/top2)-1)*100:.1f}% premium for top pathway over second")
print(f"   Entropy captures this as log2 difference: {np.log2(top1/top2):.3f} bits")

# D3: Distributional consistency
# Add a weak pathway (|z| = 0.5) and see how each score changes
vals_perturbed = np.append(vals, 0.5)
s_p = s_score(vals_perturbed, threshold=2.0)
h_p = entropy_score(vals_perturbed)
g_p = gini_score(vals_perturbed)
r_p = ratio_score(vals_perturbed)

print(f"\nD3 (distributional consistency): adding weak pathway (|z|=0.5)")
print(f"   S-score:  {s:.4f} -> {s_p:.4f}  (change: {s_p-s:+.4f})")
print(f"   Entropy:  {h:.4f} -> {h_p:.4f}  (change: {h_p-h:+.4f})")
print(f"   Gini:     {g:.4f} -> {g_p:.4f}  (change: {g_p-g:+.4f})")
print(f"   Ratio:    {r:.4f} -> {r_p:.4f}  (change: {r_p-r:+.4f})")
print(f"   Ratio is {'stable' if abs(r_p-r) < 0.01 else 'unstable'} (unaffected by weak pathway addition)")
print(f"   Entropy {'increases' if h_p > h else 'decreases'} (becomes {'less' if h_p > h else 'more'} specific)")

# D4: Monotonicity under weak off-target addition
# If we strengthen a near-zero pathway slightly, does specificity decrease?
vals_d4 = vals.copy()
# Find the weakest pathway and strengthen it slightly
weakest_idx = np.argmin(vals_d4)
vals_d4[weakest_idx] += 1.0  # add signal to weakest pathway

s_d4 = s_score(vals_d4, threshold=2.0)
h_d4 = entropy_score(vals_d4)
g_d4 = gini_score(vals_d4)
r_d4 = ratio_score(vals_d4)

print(f"\nD4 (monotonicity): strengthening weakest pathway by |z|=1.0")
print(f"   S-score:  {s:.4f} -> {s_d4:.4f}  {'✓ decreases (correct)' if s_d4 >= s else '✗ increases (violation)'}")
print(f"   Entropy:  {h:.4f} -> {h_d4:.4f}  {'✓ increases (correct)' if h_d4 >= h else '✗ decreases (violation)'}")
print(f"   Gini:     {g:.4f} -> {g_d4:.4f}  {'✓ decreases (correct)' if g_d4 <= g else '✗ increases (violation)'}")
print(f"   Ratio:    {r:.4f} -> {r_d4:.4f}  {'✓ unchanged (correct)' if abs(r_d4-r) < 0.01 else '✗ changes (violation)'}")

# ── 5. Sensitivity analysis ───────────────────────────────────────────────────
# How do scores change as we vary the threshold for D1?
thresholds = np.arange(0.5, 5.0, 0.5)
s_scores = [s_score(vals, t) for t in thresholds]

# How do scores change as we add progressively stronger off-target pathways?
extra_signals = np.arange(0, 5.0, 0.25)
entropies  = [entropy_score(np.append(vals, e)) for e in extra_signals]
ginis      = [gini_score(np.append(vals, e)) for e in extra_signals]
ratios     = [ratio_score(np.append(vals, e)) for e in extra_signals]
s_scores2  = [s_score(np.append(vals, e), threshold=2.0) for e in extra_signals]

# ── 6. Plot ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 3, figure=fig)

# Panel A: pathway z-scores
ax1 = fig.add_subplot(gs[0, :2])
colors = ['tomato' if v > 0 else 'steelblue' for v in progeny.sort_values().values]
ax1.barh(range(len(progeny)), progeny.sort_values().values, color=colors)
ax1.set_yticks(range(len(progeny)))
ax1.set_yticklabels(progeny.sort_values().index, fontsize=9)
ax1.axvline(0, color='black', linewidth=0.8)
ax1.axvline(2, color='grey', linestyle=':', alpha=0.7)
ax1.axvline(-2, color='grey', linestyle=':', alpha=0.7)
ax1.set_xlabel('PROGENy z-score')
ax1.set_title('A. Pathway activity profile (imatinib vs control, K562)', fontweight='bold')

# Panel B: scores table
ax2 = fig.add_subplot(gs[0, 2])
ax2.axis('off')
table_data = [
    ['Metric', 'Value', 'Interpretation'],
    ['S-score', f'{s:.3f}', f'{int(s*14)}/14 pathways |z|>2'],
    ['Entropy', f'{h:.3f}', 'bits (lower=specific)'],
    ['Gini', f'{g:.3f}', 'higher=specific'],
    ['Ratio', f'{r:.3f}', 'top/second pathway'],
]
table = ax2.table(cellText=table_data[1:], colLabels=table_data[0],
                  loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.8)
ax2.set_title('B. Transcriptional\nspecificity scores', fontweight='bold')

# Panel C: D3 — response to adding weak pathways
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(extra_signals, entropies, 'o-', color='steelblue', label='Entropy', markersize=4)
ax3_twin = ax3.twinx()
ax3_twin.plot(extra_signals, ginis, 's--', color='tomato', label='Gini', markersize=4)
ax3.set_xlabel('Added pathway |z|')
ax3.set_ylabel('Entropy (bits)', color='steelblue')
ax3_twin.set_ylabel('Gini', color='tomato')
ax3.set_title('C. D3: Adding off-target pathway', fontweight='bold')
ax3.legend(loc='upper left', fontsize=8)
ax3_twin.legend(loc='upper right', fontsize=8)

# Panel D: ratio response
ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(extra_signals, ratios, 'o-', color='purple', markersize=4)
ax4.axhline(r, color='grey', linestyle=':', label='Original ratio')
ax4.set_xlabel('Added pathway |z|')
ax4.set_ylabel('Ratio (top/second)')
ax4.set_title('D. D3: Ratio stability', fontweight='bold')
ax4.legend(fontsize=8)

# Panel E: S-score threshold sensitivity
ax5 = fig.add_subplot(gs[1, 2])
ax5.plot(thresholds, s_scores, 'o-', color='green', markersize=6)
ax5.set_xlabel('|z| threshold')
ax5.set_ylabel('S-score (fraction perturbed)')
ax5.set_title('E. D1: Threshold sensitivity', fontweight='bold')
ax5.grid(True, alpha=0.3)

plt.suptitle('Transcriptional Specificity Analysis — Imatinib (K562)\nAnalog of binding selectivity desiderata',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('analysis/transcriptional_specificity/specificity_analysis.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nOutput: analysis/transcriptional_specificity/specificity_analysis.png")
