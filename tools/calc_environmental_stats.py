#!/usr/bin/env python3
"""
calc_environmental_stats.py - Reproduces Environmental Statistics for Report (Page 6)

Calculates:
1. Oceanographic metadata distribution (depth layers, temperature range)
2. Dominant PETase homolog clusters and their environmental niches
3. Total PETase abundance across vertical depth layers + Kruskal-Wallis test
4. Spearman rank correlations between PETase abundance and environmental variables
"""

import pandas as pd
import numpy as np
from scipy import stats

def main():
    # 1. Load metadata and abundance tables
    meta = pd.read_csv('data/metadata.csv', index_col=0)
    ab = pd.read_csv('data/ab_hits.tsv', sep='\t', index_col=0)

    # Align stations present in both datasets
    common = [c for c in meta.index if c in ab.columns]
    meta = meta.loc[common].copy()
    ab = ab[common]

    print("============================================================")
    print("1. TARA OCEANS STATIONS & ENVIRONMENTAL OVERVIEW")
    print("============================================================")
    print(f"Total sampled stations: {len(meta)}")
    print(f"Water temperature range: {meta.temperature.min():.1f} °C to {meta.temperature.max():.1f} °C")

    # Define vertical depth layers
    meta['depth_layer'] = pd.cut(
        meta['depth'],
        bins=[-1, 10, 150, 2000],
        labels=['Surface (0-10m)', 'DCM (10-150m)', 'Mesopelagic (>150m)']
    )
    print("\nStation counts by depth layer:")
    for layer, count in meta['depth_layer'].value_counts()[['Surface (0-10m)', 'DCM (10-150m)', 'Mesopelagic (>150m)']].items():
        pct = count / len(meta) * 100
        print(f"  - {layer:20s}: {count:2d} stations ({pct:.1f}%)")

    # Total PETase abundance
    meta['total_petase'] = ab.sum(axis=0)
    detected = (meta['total_petase'] > 0).sum()
    print(f"\nPETase detected in {detected} of {len(meta)} stations.")

    print("\n============================================================")
    print("2. DOMINANT PETASE HOMOLOG CLUSTERS (NICHE PARTITIONING)")
    print("============================================================")
    meta['dominant'] = ab.idxmax(axis=0)
    top4 = meta['dominant'].value_counts().nlargest(4)
    total_top4 = top4.sum()
    print(f"Top 4 clusters dominate {total_top4} of {len(meta)} stations:\n")

    for cl, count in top4.items():
        sub = meta[meta['dominant'] == cl]
        cid = cl[-9:]
        t_mean, t_std = sub.temperature.mean(), sub.temperature.std()
        d_mean, d_std = sub.depth.mean(), sub.depth.std()
        lat_min, lat_max = sub['geographic location (latitude)'].min(), sub['geographic location (latitude)'].max()
        print(f"Cluster {cid} (n={count:2d} stations):")
        print(f"  - Temperature: {t_mean:4.1f} ± {t_std:.1f} °C")
        print(f"  - Depth:       {d_mean:5.1f} ± {d_std:.1f} m")
        print(f"  - Latitudes:   {lat_min:5.1f}° to {lat_max:5.1f}°")

    print("\n============================================================")
    print("3. VERTICAL STRATIFICATION ACROSS OCEAN LAYERS")
    print("============================================================")
    layer_order = ['Surface (0-10m)', 'DCM (10-150m)', 'Mesopelagic (>150m)']
    for layer in layer_order:
        vals = meta[meta['depth_layer'] == layer]['total_petase'].dropna()
        print(f"{layer:20s}: Mean = {vals.mean():5.1f} ± {vals.std():.1f}, Median = {vals.median():5.1f} (n={len(vals)})")

    groups = [meta[meta['depth_layer'] == l]['total_petase'].dropna() for l in layer_order]
    H, p_kw = stats.kruskal(*groups)
    print(f"\nNon-parametric Kruskal-Wallis test: H = {H:.2f}, p = {p_kw:.2e}")

    print("\n============================================================")
    print("4. SPEARMAN RANK CORRELATIONS (TOTAL ABUNDANCE VS ENVIRONMENT)")
    print("============================================================")
    variables = [
        ('temperature', 'Water Temperature (°C)'),
        ('depth', 'Hydrostatic Depth (m)'),
        ('salinity', 'Salinity (PSU)'),
        ('oxygen sensor', 'Dissolved Oxygen (μmol/kg)'),
        ('chlorophyll sensor', 'Chlorophyll Sensor (mg/m³)')
    ]
    for col, desc in variables:
        sub = meta.dropna(subset=[col, 'total_petase'])
        rho, p = stats.spearmanr(sub['total_petase'], sub[col])
        print(f"{desc:30s}: Spearman rho = {rho:+.3f}, p = {p:.2e} (n={len(sub)})")
    print("============================================================")

if __name__ == '__main__':
    main()

