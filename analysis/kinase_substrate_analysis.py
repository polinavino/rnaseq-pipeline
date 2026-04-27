import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, pearsonr
import urllib.request
import json
import time

# ── 1. Load imatinib binding profile ─────────────────────────────────────────
klaeger = pd.read_csv(os.environ.get("KLAEGER_PATH", "/Users/polina/Documents/BioInfStuff/selectivity/klaeger_matrix.csv"), index_col=0)
imatinib = klaeger.loc['Imatinib'].dropna().sort_values(ascending=False)
print(f"Imatinib binding profile: {len(imatinib)} kinases")
print(f"Top 10 targets:\n{imatinib.head(10)}\n")

# ── 2. Load DESeq2 results ────────────────────────────────────────────────────
deseq2 = pd.read_csv('analysis/deseq2_results.csv').set_index('gene_id')
deseq2 = deseq2.dropna(subset=['log2FoldChange', 'padj'])
print(f"DESeq2 results: {len(deseq2)} genes")

# ── 3. Fetch kinase-substrate relationships from Harmonizome ─────────────────
def get_substrates(kinase):
    url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{kinase}/PhosphoSitePlus+Substrates+of+Kinases"
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.loads(response.read())
            substrates = [a['gene']['symbol'] for a in data.get('associations', [])]
            return substrates
    except:
        return []

# Fetch substrates for top 50 imatinib targets
print("Fetching kinase-substrate relationships from Harmonizome...")
results = []
top_kinases = imatinib.head(50).index.tolist()

for i, kinase in enumerate(top_kinases):
    substrates = get_substrates(kinase)
    time.sleep(0.3)  # rate limit
    
    if len(substrates) == 0:
        continue
    
    # Get expression changes for substrates
    substrate_lfc = []
    for sub in substrates:
        if sub in deseq2.index:
            substrate_lfc.append(deseq2.loc[sub, 'log2FoldChange'])
    
    if len(substrate_lfc) < 3:
        continue
    
    results.append({
        'kinase': kinase,
        'pKd': imatinib[kinase],
        'n_substrates': len(substrates),
        'n_substrates_detected': len(substrate_lfc),
        'mean_substrate_lfc': np.mean(substrate_lfc),
        'mean_abs_substrate_lfc': np.mean(np.abs(substrate_lfc)),
        'median_substrate_lfc': np.median(substrate_lfc),
    })
    
    if (i + 1) % 10 == 0:
        print(f"  Processed {i+1}/{len(top_kinases)} kinases...")

df = pd.DataFrame(results)
print(f"\nKinases with sufficient substrate data: {len(df)}")
print(df.to_string(index=False))

# ── 4. Correlate binding affinity with substrate response ─────────────────────
r_mean, p_mean = spearmanr(df['pKd'], df['mean_substrate_lfc'])
r_abs, p_abs   = spearmanr(df['pKd'], df['mean_abs_substrate_lfc'])

print(f"\nSpearman correlation (pKd vs mean substrate log2FC): r={r_mean:.3f}, p={p_mean:.4f}")
print(f"Spearman correlation (pKd vs mean |substrate log2FC|): r={r_abs:.3f}, p={p_abs:.4f}")

# ── 5. Key targets ────────────────────────────────────────────────────────────
print("\nKey imatinib targets:")
for target in ['ABL1', 'ABL2', 'KIT', 'PDGFRB', 'PDGFRA', 'DDR1', 'DDR2']:
    row = df[df['kinase'] == target]
    if len(row) > 0:
        r = row.iloc[0]
        print(f"  {target:8s}: pKd={r['pKd']:.2f}, n_substrates={r['n_substrates_detected']:3.0f}, "
              f"mean_lfc={r['mean_substrate_lfc']:+.3f}, mean_|lfc|={r['mean_abs_substrate_lfc']:.3f}")

# ── 6. Save results ───────────────────────────────────────────────────────────
df.to_csv('analysis/kinase_substrate_results.csv', index=False)

# ── 7. Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: pKd vs mean substrate log2FC
axes[0].scatter(df['pKd'], df['mean_substrate_lfc'], 
                s=df['n_substrates_detected']*3, alpha=0.7, color='steelblue')
axes[0].axhline(0, color='grey', linestyle='--', alpha=0.5)

for _, row in df[df['kinase'].isin(['ABL1','ABL2','KIT','PDGFRB','DDR1'])].iterrows():
    axes[0].annotate(row['kinase'], (row['pKd'], row['mean_substrate_lfc']),
                     fontsize=8, xytext=(5,5), textcoords='offset points')

axes[0].set_xlabel('Binding affinity (pKd)')
axes[0].set_ylabel('Mean log2FC of substrates')
axes[0].set_title(f'A. Binding affinity vs substrate transcriptional response\nSpearman r={r_mean:.3f}, p={p_mean:.4f}',
                  fontweight='bold')

# Panel B: pKd vs mean |substrate log2FC|
axes[1].scatter(df['pKd'], df['mean_abs_substrate_lfc'],
                s=df['n_substrates_detected']*3, alpha=0.7, color='tomato')

for _, row in df[df['kinase'].isin(['ABL1','ABL2','KIT','PDGFRB','DDR1'])].iterrows():
    axes[1].annotate(row['kinase'], (row['pKd'], row['mean_abs_substrate_lfc']),
                     fontsize=8, xytext=(5,5), textcoords='offset points')

axes[1].set_xlabel('Binding affinity (pKd)')
axes[1].set_ylabel('Mean |log2FC| of substrates')
axes[1].set_title(f'B. Binding affinity vs substrate response magnitude\nSpearman r={r_abs:.3f}, p={p_abs:.4f}',
                  fontweight='bold')

plt.suptitle('Kinase-Substrate Analysis: Does binding affinity predict substrate transcriptional response?\nImatinib (K562), PhosphoSitePlus substrates, DESeq2 results',
             fontsize=10, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('analysis/kinase_substrate_analysis.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nOutput: analysis/kinase_substrate_analysis.png")
print("Output: analysis/kinase_substrate_results.csv")
