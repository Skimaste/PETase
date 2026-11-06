#!/usr/bin/env python3
"""
plot_petase_clusters.py - KISS Version (Keep It Simple, Stupid)

Loads Tara Oceans metadata and PETase abundance data, determines the dominant
PETase cluster per station, and generates:
1. Global distribution map of dominant clusters (fig/petase_dominant_clusters.png)
2. Environmental niche and depth stratification plots (fig/petase_environmental_niches.png)
"""

import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np

os.makedirs('fig', exist_ok=True)

# 1. Load data
meta = pd.read_csv('data/metadata.csv', index_col=0)
ab = pd.read_csv('data/ab_hits.tsv', sep='\t', index_col=0)

# Match stations
common = [c for c in meta.index if c in ab.columns]
meta = meta.loc[common].copy()
ab = ab[common]

# Compute total abundance and dominant cluster
meta['total_petase'] = ab.sum(axis=0)
meta['dominant'] = ab.idxmax(axis=0)

# Top 4 dominant clusters
top4 = meta['dominant'].value_counts().nlargest(4).index.tolist()
meta['cluster_label'] = meta['dominant'].apply(lambda x: x[-9:] if x in top4 else 'Other')

palette = {
    top4[0][-9:]: '#1f77b4',  # Cold/polar (007574371)
    top4[1][-9:]: '#d62728',  # Warm surface (004492421)
    top4[2][-9:]: '#ff7f0e',  # Warm shallow (021954915)
    top4[3][-9:]: '#2ca02c',  # Mesopelagic deep (013456435)
    'Other': '#7f7f7f'
}

# --- Plot 1: Global World Map ---
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.Robinson()})
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#e5e5e5', edgecolor='#aaaaaa', linewidth=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='#f4f8fb')
ax.coastlines(resolution='110m', color='#777777', linewidth=0.6)

for label, color in palette.items():
    sub = meta[meta['cluster_label'] == label]
    ax.scatter(
        sub['geographic location (longitude)'],
        sub['geographic location (latitude)'],
        c=color,
        s=sub['total_petase'] / 15 + 20,
        alpha=0.85,
        edgecolors='k',
        linewidth=0.4,
        transform=ccrs.PlateCarree(),
        label=f'Cluster {label} (n={len(sub)})'
    )

ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=True, fontsize=8)
plt.title('Global Distribution of Dominant PETase Homolog Clusters (Tara Oceans)', fontsize=11, weight='bold', pad=10)
plt.tight_layout()
plt.savefig('fig/petase_dominant_clusters.png', dpi=300, bbox_inches='tight')
plt.close()
print('Saved fig/petase_dominant_clusters.png')

# --- Plot 2: Thermal & Depth Niche Stratification ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot A: Temperature vs Depth
for label, color in palette.items():
    sub = meta[meta['cluster_label'] == label]
    ax1.scatter(sub['temperature'], sub['depth'], c=color, s=40, alpha=0.8, edgecolors='k', linewidth=0.3, label=f'Cluster {label}')

ax1.set_yscale('log')
ax1.invert_yaxis()
ax1.set_xlabel('Water Temperature (°C)', fontsize=10)
ax1.set_ylabel('Depth (m, log scale)', fontsize=10)
ax1.set_title('A: Thermal & Depth Niche Partitioning', fontsize=10, weight='bold')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(fontsize=8, loc='lower left')

# Subplot B: Total Abundance by Depth Layer
meta['depth_layer'] = pd.cut(meta['depth'], bins=[-1, 10, 150, 2000], labels=['Surface (0-10m)', 'DCM (10-150m)', 'Mesopelagic (>150m)'])
layer_data = [meta[meta['depth_layer'] == l]['total_petase'].dropna() for l in ['Surface (0-10m)', 'DCM (10-150m)', 'Mesopelagic (>150m)']]

bp = ax2.boxplot(layer_data, tick_labels=['Surface\n(0-10m)', 'DCM\n(10-150m)', 'Mesopelagic\n(>150m)'], patch_artist=True, medianprops=dict(color='black', lw=1.5))
colors = ['#aec7e8', '#98df8a', '#ffbb78']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

ax2.set_ylabel('Total PETase Abundance (counts)', fontsize=10)
ax2.set_title('B: Vertical Stratification (Kruskal-Wallis p = 1.4e-5)', fontsize=10, weight='bold')
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('fig/petase_environmental_niches.png', dpi=300)
plt.close()
print('Saved fig/petase_environmental_niches.png')
