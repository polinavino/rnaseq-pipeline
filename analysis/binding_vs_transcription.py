import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ── 1. Load Klaeger binding data for imatinib ─────────────────────────────────
klaeger = pd.read_csv('/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv', index_col=0)

imatinib_binding = klaeger.loc['Imatinib'].dropna().sort_values(ascending=False)
print(f"Imatinib binds {len(imatinib_binding)} kinases")
print("\nTop 20 imatinib targets by pKd:")
print(imatinib_binding.head(20))

# ── 2. Load DESeq2 results ────────────────────────────────────────────────────
deseq2 = pd.read_csv('analysis/deseq2_results.csv')
deseq2 = deseq2.dropna(subset=['padj', 'log2FoldChange'])

# ── 3. Map kinases to gene expression changes ─────────────────────────────────
# Check which of imatinib's kinase targets are differentially expressed
kinase_expr = []
for kinase in imatinib_binding.index:
    match = deseq2[deseq2['gene_id'] == kinase]
    if len(match) > 0:
        row = match.iloc[0]
        kinase_expr.append({
            'kinase': kinase,
            'pKd': imatinib_binding[kinase],
            'log2FC': row['log2FoldChange'],
            'padj': row['padj'],
            'significant': row['padj'] < 0.05
        })

kinase_df = pd.DataFrame(kinase_expr)
print(f"\n{len(kinase_df)} of {len(imatinib_binding)} imatinib targets found in RNA-seq data")
print(f"Significantly DE: {kinase_df['significant'].sum()}")

# ── 4. Correlation: binding affinity vs expression change ─────────────────────
r, p = stats.pearsonr(kinase_df['pKd'], kinase_df['log2FC'])
print(f"\nCorrelation (pKd vs log2FC): r = {r:.3f}, p = {p:.4f}")

# ── 5. Plot binding affinity vs transcriptional change ───────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Scatter plot
colors = kinase_df['significant'].map({True: 'tomato', False: 'steelblue'})
axes[0].scatter(kinase_df['pKd'], kinase_df['log2FC'], 
                c=colors, alpha=0.7, s=50)

# Label top targets
top_targets = ['ABL1', 'ABL2', 'KIT', 'PDGFRA', 'PDGFRB', 'BCR']
for _, row in kinase_df[kinase_df['kinase'].isin(top_targets)].iterrows():
    axes[0].annotate(row['kinase'], 
                     (row['pKd'], row['log2FC']),
                     fontsize=8, ha='left',
                     xytext=(5, 5), textcoords='offset points')

axes[0].axhline(0, color='grey', linestyle='--', alpha=0.5)
axes[0].axhline(1, color='grey', linestyle=':', alpha=0.5)
axes[0].axhline(-1, color='grey', linestyle=':', alpha=0.5)
axes[0].set_xlabel('Binding affinity (pKd)')
axes[0].set_ylabel('log2 Fold Change (imatinib vs control)')
axes[0].set_title(f'Binding affinity vs transcriptional change\nr = {r:.3f}, p = {p:.4f}')

sig_patch = mpatches.Patch(color='tomato', label='padj < 0.05')
ns_patch = mpatches.Patch(color='steelblue', label='not significant')
axes[0].legend(handles=[sig_patch, ns_patch])

# Bar plot of top 15 targets by binding affinity
top15 = kinase_df.nlargest(15, 'pKd')
bar_colors = top15['log2FC'].apply(lambda x: 'tomato' if x > 0 else 'steelblue')
axes[1].barh(range(len(top15)), top15['log2FC'], color=bar_colors, alpha=0.8)
axes[1].set_yticks(range(len(top15)))
axes[1].set_yticklabels(top15['kinase'], fontsize=9)
axes[1].axvline(0, color='black', linewidth=0.8)
axes[1].axvline(1, color='grey', linestyle=':', alpha=0.5)
axes[1].axvline(-1, color='grey', linestyle=':', alpha=0.5)
axes[1].set_xlabel('log2 Fold Change')
axes[1].set_title('Top 15 imatinib targets by binding affinity\n(transcriptional change in K562)')

for i, (_, row) in enumerate(top15.iterrows()):
    axes[1].text(row['log2FC'] + 0.05 if row['log2FC'] >= 0 else row['log2FC'] - 0.05,
                 i, f"pKd={row['pKd']:.1f}", 
                 va='center', ha='left' if row['log2FC'] >= 0 else 'right',
                 fontsize=7)

plt.tight_layout()
plt.savefig('analysis/binding_vs_transcription.png', dpi=150, bbox_inches='tight')
plt.close()

# ── 6. Summary ────────────────────────────────────────────────────────────────
print("\n── Top 15 imatinib targets: binding vs transcription ──")
print(top15[['kinase', 'pKd', 'log2FC', 'padj', 'significant']].to_string(index=False))

print("\n── Key targets ──")
for target in ['ABL1', 'KIT', 'PDGFRA', 'PDGFRB']:
    row = kinase_df[kinase_df['kinase'] == target]
    if len(row) > 0:
        r = row.iloc[0]
        print(f"{target:8s}: pKd={r['pKd']:.2f}, log2FC={r['log2FC']:.3f}, padj={r['padj']:.2e}, sig={r['significant']}")

print("\nOutput: analysis/binding_vs_transcription.png")
