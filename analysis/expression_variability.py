import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.stats import spearmanr, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# ── 1. Load count data ────────────────────────────────────────────────────────
counts = pd.read_csv('results/salmon/salmon.merged.gene_counts.tsv', sep='\t')
counts = counts.set_index('gene_id').drop(columns=['gene_name'])

print("Count matrix shape:", counts.shape)
print("Columns:", counts.columns.tolist())

# Separate conditions
ctrl_cols = [c for c in counts.columns if 'control' in c]
ima_cols  = [c for c in counts.columns if 'imatinib' in c]

print(f"\nControl replicates: {ctrl_cols}")
print(f"Imatinib replicates: {ima_cols}")

# ── 2. Filter low-count genes ─────────────────────────────────────────────────
# Keep genes with mean count >= 10 in at least one condition
ctrl_mean = counts[ctrl_cols].mean(axis=1)
ima_mean  = counts[ima_cols].mean(axis=1)
keep = (ctrl_mean >= 10) | (ima_mean >= 10)
counts = counts[keep]
ctrl_mean = ctrl_mean[keep]
ima_mean  = ima_mean[keep]
print(f"\nGenes after filtering: {len(counts)}")

# ── 3. Compute coefficient of variation ───────────────────────────────────────
# CV = std / mean — normalized measure of variability
ctrl_cv = counts[ctrl_cols].std(axis=1) / (counts[ctrl_cols].mean(axis=1) + 1)
ima_cv  = counts[ima_cols].std(axis=1)  / (counts[ima_cols].mean(axis=1) + 1)

# Change in variability after imatinib
delta_cv = ima_cv - ctrl_cv

# ── 4. Load DESeq2 results ────────────────────────────────────────────────────
deseq2 = pd.read_csv('analysis/deseq2_results.csv')
deseq2 = deseq2.set_index('gene_id').dropna(subset=['padj', 'log2FoldChange'])

# Merge
df = pd.DataFrame({
    'ctrl_mean':  ctrl_mean,
    'ima_mean':   ima_mean,
    'ctrl_cv':    ctrl_cv,
    'ima_cv':     ima_cv,
    'delta_cv':   delta_cv,
}).join(deseq2[['log2FoldChange', 'padj', 'stat']], how='inner')

df['significant'] = df['padj'] < 0.05
df['abs_lfc'] = df['log2FoldChange'].abs()

print(f"Genes in merged dataset: {len(df)}")

# ── 5. Question 1: Does control CV predict DE magnitude? ─────────────────────
# Kaern's model: high-variability genes are candidates for noise-driven selection
r, p = spearmanr(df['ctrl_cv'], df['abs_lfc'])
print(f"\nQ1: Spearman correlation (control CV vs |log2FC|): r={r:.3f}, p={p:.4f}")

# Compare CV of significant vs non-significant genes
sig_cv   = df.loc[df['significant'], 'ctrl_cv']
nonsig_cv = df.loc[~df['significant'], 'ctrl_cv']
u_stat, u_p = mannwhitneyu(sig_cv, nonsig_cv, alternative='two-sided')
print(f"Median CV — significant DEGs: {sig_cv.median():.4f}")
print(f"Median CV — non-significant:  {nonsig_cv.median():.4f}")
print(f"Mann-Whitney U test: p={u_p:.4e}")

# ── 6. Question 2: Does imatinib reduce expression variability? ───────────────
# Prediction: imatinib selects for uniform cells, reducing population-level CV
mean_ctrl_cv = ctrl_cv.mean()
mean_ima_cv  = ima_cv.mean()
print(f"\nQ2: Mean CV — control:  {mean_ctrl_cv:.4f}")
print(f"    Mean CV — imatinib: {mean_ima_cv:.4f}")
print(f"    Change: {mean_ima_cv - mean_ctrl_cv:+.4f} ({'reduction' if mean_ima_cv < mean_ctrl_cv else 'increase'})")

paired_t, paired_p = stats.ttest_rel(ctrl_cv, ima_cv)
print(f"    Paired t-test: t={paired_t:.3f}, p={paired_p:.4e}")

# ── 7. Question 3: Which high-CV control genes are most DE? ──────────────────
# Top 10% by control CV
high_cv_threshold = df['ctrl_cv'].quantile(0.90)
high_cv_genes = df[df['ctrl_cv'] >= high_cv_threshold]
sig_high_cv = high_cv_genes[high_cv_genes['significant']]

print(f"\nQ3: High-CV genes (top 10%): {len(high_cv_genes)}")
print(f"    Significant among high-CV: {len(sig_high_cv)} ({100*len(sig_high_cv)/len(high_cv_genes):.1f}%)")
print(f"    Significant overall: {df['significant'].mean()*100:.1f}%")

print("\nTop 10 high-CV significant genes by |log2FC|:")
top = sig_high_cv.nlargest(10, 'abs_lfc')[['ctrl_cv', 'log2FoldChange', 'padj']]
print(top.round(4))

# ── 8. Plot ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 10))
gs = gridspec.GridSpec(2, 3, figure=fig)

# Panel A: control CV vs |log2FC|
ax1 = fig.add_subplot(gs[0, :2])
colors = df['significant'].map({True: 'tomato', False: 'steelblue'})
ax1.scatter(df['ctrl_cv'], df['abs_lfc'], c=colors, alpha=0.2, s=5)
ax1.set_xlabel('Coefficient of variation (control replicates)')
ax1.set_ylabel('|log2 Fold Change| (imatinib vs control)')
ax1.set_title(f'A. Control expression variability vs transcriptional response\nSpearman r={r:.3f}, p={p:.4f}',
              fontweight='bold')
from matplotlib.patches import Patch
ax1.legend(handles=[Patch(color='tomato', label='padj<0.05'),
                    Patch(color='steelblue', label='not significant')],
           fontsize=9)

# Panel B: CV distribution
ax2 = fig.add_subplot(gs[0, 2])
ax2.hist(ctrl_cv, bins=50, alpha=0.6, color='steelblue', label='Control', density=True)
ax2.hist(ima_cv,  bins=50, alpha=0.6, color='tomato', label='Imatinib', density=True)
ax2.axvline(mean_ctrl_cv, color='steelblue', linestyle='--', linewidth=1.5)
ax2.axvline(mean_ima_cv,  color='tomato', linestyle='--', linewidth=1.5)
ax2.set_xlabel('Coefficient of variation')
ax2.set_ylabel('Density')
ax2.set_title(f'B. CV distribution shift\nMean: ctrl={mean_ctrl_cv:.3f}, ima={mean_ima_cv:.3f}',
              fontweight='bold')
ax2.legend(fontsize=9)
ax2.set_xlim(0, 2)

# Panel C: delta CV distribution
ax3 = fig.add_subplot(gs[1, 0])
ax3.hist(delta_cv, bins=50, color='purple', alpha=0.7)
ax3.axvline(0, color='black', linestyle='--')
ax3.axvline(delta_cv.mean(), color='red', linestyle='-', linewidth=1.5,
            label=f'Mean={delta_cv.mean():.3f}')
ax3.set_xlabel('ΔCVV (imatinib − control)')
ax3.set_ylabel('Number of genes')
ax3.set_title('C. Change in variability per gene', fontweight='bold')
ax3.legend(fontsize=9)

# Panel D: CV of significant vs non-significant
ax4 = fig.add_subplot(gs[1, 1])
ax4.boxplot([nonsig_cv.values, sig_cv.values],
            labels=['Not significant', 'Significant DEGs'],
            patch_artist=True,
            boxprops=dict(facecolor='lightblue'),
            medianprops=dict(color='red', linewidth=2))
ax4.set_ylabel('Control CV')
ax4.set_title(f'D. Control CV by DE status\nMann-Whitney p={u_p:.2e}', fontweight='bold')

# Panel E: High-CV genes enrichment
ax5 = fig.add_subplot(gs[1, 2])
categories = ['All genes', 'High-CV genes\n(top 10%)']
sig_fractions = [df['significant'].mean() * 100,
                 len(sig_high_cv) / len(high_cv_genes) * 100]
bars = ax5.bar(categories, sig_fractions, color=['steelblue', 'tomato'], alpha=0.8)
ax5.set_ylabel('% significantly DE (padj < 0.05)')
ax5.set_title('E. DE enrichment in high-CV genes', fontweight='bold')
for bar, val in zip(bars, sig_fractions):
    ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{val:.1f}%', ha='center', fontweight='bold')
ax5.set_ylim(0, max(sig_fractions) * 1.2)

plt.suptitle('Expression Variability Analysis — Imatinib Response (K562)\nMotivated by Kaern et al.: noise-driven drug resistance framework',
             fontsize=11, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('analysis/expression_variability.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nOutput: analysis/expression_variability.png")
