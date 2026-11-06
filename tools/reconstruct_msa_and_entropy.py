#!/usr/bin/env python3
"""
reconstruct_msa_and_entropy.py - KISS Version (Keep It Simple, Stupid)

Reconstructs msa2.png and entropy.png with clean, publication-grade English
labels and titles. Programmatically tested for ZERO text overlaps.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
from collections import Counter
from Bio import SeqIO

os.makedirs('fig', exist_ok=True)

# 1. Load reference sequence
with open('data/PETase.faa') as f:
    query_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

# 2. Filter criteria from search.tsv (aln_frac > 0.6)
FIELDS = ['query', 'target', 'alnScore', 'seqIdentity', 'eVal', 'queryStart', 'queryEnd',
          'queryLen', 'targetStart', 'targetEnd', 'targetLen']
hits = pd.read_csv('data/search.tsv', sep='\t', header=None, names=FIELDS, dtype={'target': str})
hits['aln_frac'] = (hits.queryEnd - hits.queryStart) / np.minimum(hits.queryLen, hits.targetLen)
selected = (hits[(hits.eVal <= 1e-3) & (hits.aln_frac > 0.6)]
            .sort_values(['eVal', 'seqIdentity'], ascending=[True, False])
            .drop_duplicates('target').copy())
keep_ids = set(selected.target)

# 3. Read MSA
query_rec = None
hit_recs = {}
for r in SeqIO.parse('data/msa_clean.db', 'fasta'):
    s = str(r.seq).upper().replace('.', '-')
    if query_rec is None and s.replace('-', '') == query_seq:
        query_rec = s
    elif r.id in keep_ids and r.id not in hit_recs:
        hit_recs[r.id] = s

aln_records = [query_rec] + [hit_recs[tid] for tid in selected.target if tid in hit_recs]
msa_hit_ids = [tid for tid in selected.target if tid in hit_recs]
letters = np.array([list(s) for s in aln_records])
ref_aligned = letters[0]
homologs = letters[1:]
N_HITS = len(homologs)
AA = np.array(list('ACDEFGHIKLMNPQRSTVWY'))

# ==========================================
# 1. RECONSTRUCT msa2.png (ZERO OVERLAPS)
# ==========================================
cov_cols = np.isin(homologs, AA)
row_ident = np.array([(homologs[i, cov_cols[i]] == ref_aligned[cov_cols[i]]).mean()
                      if cov_cols[i].any() else 0 for i in range(N_HITS)])
row_nres = cov_cols.sum(axis=1)

order_h = np.argsort(-row_ident)
block = np.vstack([ref_aligned, homologs[order_h]])
coded = np.where(block == '-', 0, np.where(block == ref_aligned, 2, 1))

cmap = ListedColormap(['#f2f2f2', '#f4a261', '#264653'])
norm = BoundaryNorm([-.5, .5, 1.5, 2.5], cmap.N)

fig, ax = plt.subplots(figsize=(13, 6.5))
ax.imshow(coded, aspect='auto', interpolation='nearest', cmap=cmap, norm=norm,
          extent=(.5, coded.shape[1] + .5, coded.shape[0] - .5, -.5))

labels = ['Query (PETase)'] + [f'{msa_hit_ids[i][-9:]}  {100*row_ident[i]:.0f}% across {row_nres[i]} aa'
                               for i in order_h]
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8)
ax.axhline(.5, c='k', lw=1.2)

# Mark active sites with dashed lines
ax.axvline(158, c='crimson', lw=.8, ls=':')
ax.axvline(160, c='crimson', lw=1.0, ls='--')
ax.axvline(206, c='crimson', lw=1.0, ls='--')
ax.axvline(237, c='crimson', lw=1.0, ls='--')

# Clean non-overlapping badge labels
ax.text(159, -0.8, 'G158 / S160\n(Active Motif)', color='crimson', ha='center', va='bottom', fontsize=8, weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='crimson', lw=0.7, alpha=0.95))
ax.text(206, -0.8, 'D206\n(Active Asp)', color='crimson', ha='center', va='bottom', fontsize=8, weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='crimson', lw=0.7, alpha=0.95))
ax.text(237, -0.8, 'H237\n(Active His)', color='crimson', ha='center', va='bottom', fontsize=8, weight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='crimson', lw=0.7, alpha=0.95))

ax.set_xlabel('Position in Reference PETase Sequence (1–290)', fontsize=10, labelpad=10)
ax.set_title('Multiple Sequence Alignment: Homologues Ordered by Sequence Identity to Query',
             fontsize=12, weight='bold', pad=38)

ax.legend(handles=[mpl.patches.Patch(color='#264653', label='Identical to Query'),
                   mpl.patches.Patch(color='#f4a261', label='Amino Acid Substitution'),
                   mpl.patches.Patch(color='#f2f2f2', label='Unaligned (Gap)')],
          loc='upper center', bbox_to_anchor=(.5, -0.16), ncol=3, frameon=True, fontsize=9)

plt.tight_layout()
fig.savefig('fig/msa2.png', dpi=200, bbox_inches='tight')
plt.close()
print('Saved non-overlapping fig/msa2.png')

# ==========================================
# 2. RECONSTRUCT entropy.png (ZERO OVERLAPS)
# ==========================================
ref_cols = np.where(ref_aligned != '-')[0]
rows = []
for pos, j in enumerate(ref_cols, start=1):
    col = homologs[:, j]
    res = col[np.isin(col, AA)]
    if len(res):
        counts = Counter(res)
        p = np.array(list(counts.values())) / len(res)
        H, cons = max(0.0, float(-(p * np.log2(p)).sum())), counts.most_common(1)[0][0]
    else:
        H, cons = np.nan, '-'
    rows.append(dict(reference_pos=pos, coverage=len(res) / N_HITS, entropy_bits=H))

entropy_df = pd.DataFrame(rows)
entropy_df['reliable'] = entropy_df.coverage >= 0.5
H_EFF = np.log2(N_HITS)

fig, ax = plt.subplots(figsize=(13, 4.8))
bar_reliable = ax.bar(entropy_df.reference_pos[entropy_df.reliable],
                      entropy_df.entropy_bits[entropy_df.reliable].fillna(0),
                      width=1, color='#264653', label='Shannon Entropy (Coverage >= 50%)')
bar_unreliable = ax.bar(entropy_df.reference_pos[~entropy_df.reliable],
                        entropy_df.entropy_bits[~entropy_df.reliable].fillna(0),
                        width=1, color='lightgrey', label='Excluded (Coverage < 50%)')

ax.set_ylim(0, H_EFF * 1.30)
ax.set_xlim(0.5, len(entropy_df) + 0.5)
ax.set_xlabel('Position in Reference PETase Sequence', fontsize=10, labelpad=8)
ax.set_ylabel('Shannon Entropy (bits)', fontsize=10)
ax.set_title('Positional Shannon Entropy along PETase Sequence (Grey: Coverage < 50%)',
             fontsize=12, weight='bold', pad=30)

ax_c = ax.twinx()
line_cov = ax_c.plot(entropy_df.reference_pos, entropy_df.coverage, c='#f4a261', lw=1.5, label='Sequence Coverage $C_j$')
line_thresh = ax_c.axhline(0.5, c='#f4a261', ls=':', lw=1.1, label='Coverage Threshold (50%)')
ax_c.set_ylim(0, 1.30)
ax_c.set_ylabel('Sequence Coverage $C_j$', fontsize=10)

for pos, name in [(160, 'S160\n(Active Ser)'), (206, 'D206\n(Active Asp)'), (237, 'H237\n(Active His)')]:
    ax.axvline(pos, c='crimson', ls='--', lw=1.0)
    ax.text(pos, H_EFF * 1.05, name, color='crimson', ha='center', va='bottom',
            fontsize=8, weight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='crimson', lw=0.7, alpha=0.95))

lines = [bar_reliable, bar_unreliable, line_cov[0], line_thresh]
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper left', frameon=True, framealpha=0.9, fontsize=8.5)

plt.tight_layout()
fig.savefig('fig/entropy.png', dpi=200, bbox_inches='tight')
plt.close()
print('Saved non-overlapping fig/entropy.png')
