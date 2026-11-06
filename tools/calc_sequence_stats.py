#!/usr/bin/env python3
"""
calc_sequence_stats.py - Reproduces Sequence & MSA Diversity Statistics for Report

Calculates:
1. MMseqs2 raw hits (E <= 1e-3) and length filtering (aln_frac > 0.60)
2. Sequence identity range of retained marine homologs
3. Positional sequence coverage across the 290 query positions
4. Positional Shannon entropy, active G-x-S-x-G motif, and catalytic triad conservation
"""

import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter

def main():
    # 1. Load reference query
    with open('data/PETase.faa') as f:
        query_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
    query_len = len(query_seq)

    print("============================================================")
    print("1. MMSEQS2 HOMOLOG SEARCH & FILTERING")
    print("============================================================")
    print(f"Reference query length: {query_len} aa (I. sakaiensis PETase)")

    FIELDS = ['query', 'target', 'alnScore', 'seqIdentity', 'eVal', 'queryStart', 'queryEnd',
              'queryLen', 'targetStart', 'targetEnd', 'targetLen']
    hits = pd.read_csv('data/search.tsv', sep='\t', header=None, names=FIELDS, dtype={'target': str})
    hits['aln_frac'] = (hits.queryEnd - hits.queryStart) / np.minimum(hits.queryLen, hits.targetLen)

    raw_hits = hits[hits.eVal <= 1e-3].drop_duplicates('target').copy()
    print(f"Raw MMseqs2 hits (E <= 1e-3): {len(raw_hits)}")

    # Filtering aln_frac > 0.60
    selected = (hits[(hits.eVal <= 1e-3) & (hits.aln_frac > 0.6)]
                .sort_values(['eVal', 'seqIdentity'], ascending=[True, False])
                .drop_duplicates('target').copy())
    print(f"High-confidence homologs (aln_frac > 0.60): {len(selected)}")
    min_id = selected.seqIdentity.min() * 100
    max_id = selected.seqIdentity.max() * 100
    print(f"Sequence identity range with query: {min_id:.1f}% to {max_id:.1f}%")

    print("\n============================================================")
    print("2. MULTIPLE SEQUENCE ALIGNMENT (MSA) & ENTROPY")
    print("============================================================")
    query_rec = None
    hit_recs = {}
    keep_ids = set(selected.target)
    for r in SeqIO.parse('data/msa_clean.db', 'fasta'):
        s = str(r.seq).upper().replace('.', '-')
        if query_rec is None and s.replace('-', '') == query_seq:
            query_rec = s
        elif r.id in keep_ids and r.id not in hit_recs:
            hit_recs[r.id] = s

    aln_records = [query_rec] + [hit_recs[tid] for tid in selected.target if tid in hit_recs]
    letters = np.array([list(s) for s in aln_records])
    ref_aligned = letters[0]
    homologs = letters[1:]
    N_HITS = len(homologs)
    print(f"Number of aligned marine homologs: {N_HITS}")
    print(f"Theoretical maximum Shannon entropy: log2({N_HITS}) = {np.log2(N_HITS):.2f} bits")

    AA = list('ACDEFGHIKLMNPQRSTVWY')
    ref_cols = np.where(ref_aligned != '-')[0]
    rows = []
    for pos, j in enumerate(ref_cols, start=1):
        col = homologs[:, j]
        res = col[np.isin(col, AA)]
        if len(res):
            counts = Counter(res)
            p = np.array(list(counts.values())) / len(res)
            H = max(0.0, float(-(p * np.log2(p)).sum()))
            consensus = counts.most_common(1)[0][0]
        else:
            H = np.nan
            consensus = '-'
        rows.append(dict(
            pos=pos,
            ref_aa=ref_aligned[j],
            consensus=consensus,
            coverage=len(res) / N_HITS,
            entropy=H
        ))

    df = pd.DataFrame(rows)
    reliable = df[df['coverage'] >= 0.5]
    print(f"Positions with >= 50% sequence coverage: {len(reliable)} of {len(df)}")
    print(f"Reliable entropy: Median = {reliable.entropy.median():.2f} bits, Max = {reliable.entropy.max():.2f} bits")

    print("\n============================================================")
    print("3. CATALYTIC MOTIF & ACTIVE TRIAD CONSERVATION")
    print("============================================================")
    triad_positions = [160, 206, 237]
    triad = df[df.pos.isin(triad_positions)]
    for _, r in triad.iterrows():
        print(f"Active site residue {r.ref_aa}{int(r.pos)}: Coverage = {r.coverage*100:4.1f}%, Entropy = {r.entropy:.2f} bits")

    # Canonical G-x-S-x-G motif (positions 158-162)
    motif_cols = [ref_cols[p-1] for p in range(158, 163)]
    cov_motif = []
    match_motif = []
    for h in homologs:
        segment = ''.join(h[motif_cols])
        if '-' not in segment and '.' not in segment:
            cov_motif.append(segment)
            if segment[0] == 'G' and segment[2] == 'S' and segment[4] == 'G':
                match_motif.append(segment)

    print(f"\nCanonical G-x-S-x-G motif (158-162):")
    print(f"  - Covering homologs: {len(cov_motif)} of {N_HITS}")
    print(f"  - Conserving exact G-x-S-x-G: {len(match_motif)} of {len(cov_motif)} ({len(match_motif)/len(cov_motif)*100:.1f}%)")
    print("============================================================")

if __name__ == '__main__':
    main()

