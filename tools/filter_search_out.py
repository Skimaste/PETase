#!/usr/bin/env python3
"""
filter_search_out.py - KISS Version (Keep It Simple, Stupid)

Filters MMseqs2 search output (search.tsv) using:
  1. E-value threshold (default <= 1e-3)
  2. Alignment coverage fraction: (queryEnd - queryStart) / min(queryLen, targetLen) > min_aln_frac (default 0.6)

Usage:
    python tools/filter_search_out.py [min_aln_frac]
"""

import sys
import pandas as pd
import numpy as np

min_frac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6

FIELDS = ['query', 'target', 'alnScore', 'seqIdentity', 'eVal',
          'queryStart', 'queryEnd', 'queryLen', 'targetStart', 'targetEnd', 'targetLen']

hits = pd.read_csv('data/search.tsv', sep='\t', header=None, names=FIELDS, dtype={'target': str})

# Compute alignment fraction over the shorter sequence
hits['aln_frac'] = (hits.queryEnd - hits.queryStart) / np.minimum(hits.queryLen, hits.targetLen)

# Filter criteria
passed = hits[(hits.eVal <= 1e-3) & (hits.aln_frac > min_frac)].copy()
passed = passed.sort_values('alnScore', ascending=False)

out_tsv = 'data/search_filtered.tsv'
passed.to_csv(out_tsv, sep='\t', index=False)

print(f"Total hits: {len(hits)}")
print(f"Filtered hits (E <= 1e-3, aln_frac > {min_frac}): {len(passed)}")
print(f"Saved to {out_tsv}")
