import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load PROGENy scores from files
ima = pd.read_csv('analysis/progeny_scores.csv').set_index('pathway')['z_score']
nil = pd.read_csv('analysis/progeny_nilotinib_scores.csv').set_index('pathway')['z_score_nilotinib']

# Align
common = ima.index.intersection(nil.index)
ima = ima.loc[common]
nil = nil.loc[common]

def transcriptional_specificity(z):
    vals = z.abs().values
    s = (vals > 2.0).sum() / len(vals)
    shifted = np.maximum(vals, 0)
    total = shifted.sum()
    p = shifted / total
    p = p[p > 0]
    h = float(-np.sum(p * np.log2(p)))
    n = len(vals)
    sorted_vals = np.sort(vals)
    indices = np.arange(1, n+1)
    g = (2 * np.sum(indices * sorted_vals)) / (n * np.sum(sorted_vals)) - (n+1)/n
    r = sorted_vals[-1] / sorted_vals[-2]
    return {'S-score': round(s,3), 'Entropy': round(h,3), 
            'Gini': round(g,3), 'Ratio': round(r,3)}

ima_scores = transcriptional_specificity(ima)
nil_scores = transcriptional_specificity(nil)

print("Imatinib transcriptional specificity:", ima_scores)
print("Nilotinib transcriptional specificity:", nil_scores)

# Plot comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pathway z-scores side by side
df = pd.DataFrame({'Imatinib': ima, 'Nilotinib': nil}).sort_values('Imatinib')
x = np.arange(len(df))
w = 0.35
axes[0].barh(x - w/2, df['Imatinib'], w, color='steelblue', alpha=0.8, label='Imatinib')
axes[0].barh(x + w/2, df['Nilotinib'], w, color='tomato', alpha=0.8, label='Nilotinib')
axes[0].set_yticks(x)
axes[0].set_yticklabels(df.index, fontsize=9)
axes[0].axvline(0, color='black', linewidth=0.8)
axes[0].legend()
axes[0].set_xlabel('PROGENy z-score')
axes[0].set_title('Pathway activity: Imatinib vs Nilotinib\n(K562 cells vs control)')

# Specificity scores
metrics = ['S-score', 'Entropy', 'Gini', 'Ratio']
ima_vals = [ima_scores[m] for m in metrics]
nil_vals = [nil_scores[m] for m in metrics]
x2 = np.arange(len(metrics))
axes[1].bar(x2 - w/2, ima_vals, w, color='steelblue', alpha=0.8, label='Imatinib')
axes[1].bar(x2 + w/2, nil_vals, w, color='tomato', alpha=0.8, label='Nilotinib')
axes[1].set_xticks(x2)
axes[1].set_xticklabels(metrics)
axes[1].legend()
axes[1].set_title('Transcriptional specificity scores\n(higher S/Entropy = less specific; higher Gini/Ratio = more specific)')
axes[1].set_ylabel('Score')

plt.tight_layout()
plt.savefig('analysis/imatinib_vs_nilotinib.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved analysis/imatinib_vs_nilotinib.png")

# Save comparison table
comp = pd.DataFrame({'imatinib': ima_scores, 'nilotinib': nil_scores})
comp.to_csv('analysis/specificity_comparison.csv')
print("Saved analysis/specificity_comparison.csv")
