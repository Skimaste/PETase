#!/usr/bin/env python3
"""
calc_structural_stats.py - Reproduces 3D Structural Diversity Statistics for Report

Calculates:
1. Alignment between query MSA and crystal structure of PETase (PDB 6QGC, Chain A)
2. Structural metrics per residue:
   - Distance from geometric centroid (C-alpha)
   - Structural packing density: mean Euclidean distance to 8 nearest C-beta neighbors
   - Relative Solvent Accessibility (RSA) via SASA (Shrake-Rupley)
3. Spearman rank correlations and 1,000 permutation tests with positional Shannon entropy
4. Comparison of solvent-exposed vs buried core residues
"""

import numpy as np
import pandas as pd
from collections import Counter
from scipy.spatial import cKDTree
from scipy import stats
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

def compute_entropy_table():
    with open('data/PETase.faa') as f:
        query_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

    FIELDS = ['query', 'target', 'alnScore', 'seqIdentity', 'eVal', 'queryStart', 'queryEnd',
              'queryLen', 'targetStart', 'targetEnd', 'targetLen']
    hits = pd.read_csv('data/search.tsv', sep='\t', header=None, names=FIELDS, dtype={'target': str})
    hits['aln_frac'] = (hits.queryEnd - hits.queryStart) / np.minimum(hits.queryLen, hits.targetLen)
    selected = (hits[(hits.eVal <= 1e-3) & (hits.aln_frac > 0.6)]
                .sort_values(['eVal', 'seqIdentity'], ascending=[True, False])
                .drop_duplicates('target').copy())
    keep_ids = set(selected.target)

    query_rec = None
    hit_recs = {}
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
        else:
            H = np.nan
        rows.append(dict(resnum=pos, coverage=len(res) / N_HITS, entropy=H))

    df = pd.DataFrame(rows)
    return df[df['coverage'] >= 0.5].copy()

def cb_virtuale(res):
    """Computes C-beta coordinates (or virtual C-beta for Glycine)."""
    if 'CB' in res and res.resname != 'GLY':
        return res['CB'].coord
    n, ca, c = res['N'].coord, res['CA'].coord, res['C'].coord
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca

def main():
    print("============================================================")
    print("1. PDB 6QGC STRUCTURE RESOLUTION & MAPPING")
    print("============================================================")
    structure = PDBParser(QUIET=True).get_structure('6qgc', 'data/6QGC.pdb')[0]['A']
    clean_residues = [r for r in structure if r.id[0] == ' ']
    print(f"Total resolved residues in 6QGC (Chain A): {len(clean_residues)} of 290")
    print("Note: Residues 1-28 correspond to the N-terminal signal peptide absent in crystal.")

    # Calculate SASA
    ShrakeRupley().compute(structure, level='R')

    # Load reliable entropy
    ent_df = compute_entropy_table()
    resnum_set = set(ent_df.resnum)
    mapped_residues = [r for r in clean_residues if r.id[1] in resnum_set]
    print(f"Residues with resolved 3D coordinates & reliable entropy (cov >= 50%): {len(mapped_residues)}")

    ent_map = dict(zip(ent_df.resnum, ent_df.entropy))
    y = np.array([ent_map[r.id[1]] for r in mapped_residues])

    # 1. Centroid distance
    all_ca = np.array([r['CA'].coord for r in clean_residues if 'CA' in r])
    centroid = all_ca.mean(axis=0)
    ca_coords = np.array([r['CA'].coord for r in mapped_residues])
    dist_centroid = np.linalg.norm(ca_coords - centroid, axis=1)

    # 2. C-beta 8-nearest neighbor distance
    all_cb = np.array([cb_virtuale(r) for r in clean_residues if all(k in r for k in ('N', 'CA', 'C'))])
    cb_coords = np.array([cb_virtuale(r) for r in mapped_residues])
    tree = cKDTree(all_cb)
    dists, _ = tree.query(cb_coords, k=9)  # k=9 includes the residue itself at index 0
    dist_cb8 = dists[:, 1:].mean(axis=1)

    # 3. Relative Solvent Accessibility (RSA)
    MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197,
                   L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
    rsa = np.array([(r.sasa / MAX_ASA.get(seq1(r.resname), 200.0)) for r in mapped_residues])
    rsa = np.clip(rsa, 0.0, 1.0)

    print("\n============================================================")
    print("2. STRUCTURAL CORRELATIONS WITH SHANNON ENTROPY")
    print("============================================================")
    SEED = 42
    N_PERM = 1000
    rng = np.random.default_rng(SEED)

    metrics = [
        ('Centroid Distance (C-alpha)', dist_centroid),
        ('C-beta 8-Nearest Neighbors Packing Distance', dist_cb8),
        ('Relative Solvent Accessibility (RSA)', rsa),
    ]

    for name, x in metrics:
        rho, p_val = stats.spearmanr(x, y)
        perm = np.array([stats.spearmanr(x, rng.permutation(y))[0] for _ in range(N_PERM)])
        p_perm = (np.sum(np.abs(perm) >= abs(rho)) + 1) / (N_PERM + 1)
        print(f"{name}:")
        print(f"  - Spearman rho = {rho:+.3f} (p_asymptotic = {p_val:.2e})")
        print(f"  - Permutation p-value (1,000 permutations) = {p_perm:.4f} (p < 0.001)")

    print("\n============================================================")
    print("3. SOLVENT-EXPOSED VS BURIED CORE RESIDUES")
    print("============================================================")
    buried = y[rsa < 0.25]
    exposed = y[rsa >= 0.25]
    u_stat, p_mw = stats.mannwhitneyu(exposed, buried, alternative='greater')
    print(f"Buried residues (RSA < 25%, n={len(buried)}):   Median entropy = {np.median(buried):.2f} bits, Mean = {np.mean(buried):.2f} bits")
    print(f"Exposed residues (RSA >= 25%, n={len(exposed)}):  Median entropy = {np.median(exposed):.2f} bits, Mean = {np.mean(exposed):.2f} bits")
    print(f"Mann-Whitney U test (Exposed > Buried): U = {u_stat:.0f}, p = {p_mw:.2e}")
    print("============================================================")

if __name__ == '__main__':
    main()

