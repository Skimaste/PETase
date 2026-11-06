#!/usr/bin/env python3
"""
plot_world_map.py - KISS Version (Keep It Simple, Stupid)

Plots oceanographic variables (temperature, salinity, depth, etc.) on a Robinson
world map for the 179 Tara Oceans stations.

Usage:
    python tools/plot_world_map.py [temperature|salinity|depth|oxygen]
Default: temperature
"""

import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np

os.makedirs('fig', exist_ok=True)

# Choose attribute to plot
attr = sys.argv[1] if len(sys.argv) > 1 else 'temperature'

# Load metadata
meta = pd.read_csv('data/metadata.csv', index_col=0)

# Attribute configurations (colormap, label, log-scale)
configs = {
    'temperature': {'col': 'temperature', 'cmap': 'coolwarm', 'label': 'Temperature (°C)', 'log': False},
    'salinity': {'col': 'salinity', 'cmap': 'viridis', 'label': 'Salinity (PSU)', 'log': False},
    'depth': {'col': 'depth', 'cmap': 'plasma_r', 'label': 'Depth (m)', 'log': True},
    'oxygen': {'col': 'oxygen sensor', 'cmap': 'cividis', 'label': 'Oxygen Sensor (μmol/kg)', 'log': False},
    'chlorophyll': {'col': 'chlorophyll sensor', 'cmap': 'YlGn', 'label': 'Chlorophyll Sensor (mg/m³)', 'log': False},
}

cfg = configs.get(attr, {'col': attr, 'cmap': 'viridis', 'label': attr, 'log': False})
val_col = cfg['col']

# Clean data
df = meta.dropna(subset=['geographic location (longitude)', 'geographic location (latitude)', val_col]).copy()
vals = df[val_col]
if cfg['log']:
    vals = np.log10(np.maximum(vals, 1.0))
    cfg['label'] += ' [log10]'

# Plot Robinson Map
fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.Robinson()})
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#eaeaea', edgecolor='#bcbcbc', linewidth=0.5)
ax.add_feature(cfeature.OCEAN, facecolor='#f7fbff')
ax.coastlines(resolution='110m', color='#888888', linewidth=0.6)

sc = ax.scatter(
    df['geographic location (longitude)'],
    df['geographic location (latitude)'],
    c=vals,
    cmap=cfg['cmap'],
    s=35,
    alpha=0.9,
    edgecolors='k',
    linewidth=0.3,
    transform=ccrs.PlateCarree()
)

cbar = plt.colorbar(sc, orientation='horizontal', pad=0.08, shrink=0.6, aspect=25)
cbar.set_label(cfg['label'], fontsize=9)
plt.title(f'Tara Oceans - Global {cfg["label"]}', fontsize=11, weight='bold', pad=10)
plt.tight_layout()

out_path = f'fig/world_map_{attr}.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()
print(f'Saved {out_path}')
